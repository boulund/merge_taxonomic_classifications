#!/usr/bin/env python3
__author__ = "Fredrik Boulund"
__date__ = "2016"
__doc__ = """Merge taxonomic read assignments produced by Kaiju, Kraken, and CLARK-S"""

from sys import argv, exit, stdout
import argparse
import sqlite3
import logging
import time
from itertools import chain, combinations, islice
from collections import namedtuple, defaultdict


def grouper(n, iterable):
    """Groups an iterable into n-sized chunks.

    :param n:  number of items per chunk.
    :param iterable:  an iterable to divide into chunks of size n.
    :return:  n-sized chunks from the iterable.
    """
    it = iter(iterable)
    while True:
       chunk = tuple(islice(it, n))
       if not chunk:
           return
       yield chunk


def parse_args():
    """Parse command line arguments.

    Returns: 
        options  argparse options namespace
        logger   logger instance
    """

    desc = "Merge outputs from Kaiju, Kraken, and CLARK-S"
    epilog = "Copyright (c) {author} {date}".format(author=__author__, date=__date__)

    parser = argparse.ArgumentParser(description=desc, epilog=epilog)

    parser.add_argument("-k", "--kaiju", dest="kaiju", metavar="KAIJU",
            help="Kaiju output file (tab).")
    parser.add_argument("-K", "--kraken", dest="kraken", metavar="KRAKEN",
            help="Kraken output file (tab).")
    parser.add_argument("-c", "--clarks", dest="clarks", metavar="CLARKS",
            help="CLARK-S output file (csv).")
    parser.add_argument("-m", "--merge-order", dest="merge_order", metavar="ORDER",
            default="clarks,kraken,kaiju",
            help="Define merge ordering, e.g. kraken,kaiju,clarks [%(default)s].")
    parser.add_argument("-d", "--dbfile", dest="dbfile", metavar="DBFILE",
            default=":memory:",
            help="Write intermediary SQLite3 database to DBFILE [%(default)s].")
    parser.add_argument("-o", "--output", dest="output", metavar="OUTFILE",
            default="merged_classifications.tab",
            help="Output filename [%(default)s].")

    parser.add_argument("--loglevel",
            choices=["DEBUG", "INFO"],
            default="DEBUG",
            help="Set log level [%(default)s].")
    parser.add_argument("--logfile", metavar="LOGFILE",
            default="",
            help="Log to LOGFILE instead of STDOUT.")

    if len(argv) < 2:
        parser.print_help()
        exit(1)

    options = parser.parse_args()

    logger = logging.getLogger(__name__)
    if options.loglevel == "DEBUG":
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)

    handlers = []
    handlers.append(logging.StreamHandler()) # Console handler
    if options.logfile:
        handlers.append(logging.FileHandler(options.logfile))

    formatter = logging.Formatter("%(asctime)s %(levelname)s: %(message)s")

    for handler in handlers:
        if options.loglevel == "DEBUG":
            handler.setLevel(logging.DEBUG)
        else:
            handler.setLevel(logging.INFO)
        handler.setFormatter(formatter)
        logger.addHandler(handler)
 
    logger.debug(" Logging started ".center(50, "="))
    return options, logger



class Merger():
    def __init__(self, dbfile):
        self.db = sqlite3.connect(dbfile)
        
        # Optimizations for fast inserts (see self.fill_table)
        self.db.execute("PRAGMA journal_mode = MEMORY")
        self.db.execute("PRAGMA synchronous = OFF")
        
        tables = ["kaiju", "clarks", "kraken", "merged"]
        drop_table = """DROP TABLE IF EXISTS {}"""
        create_table = """CREATE TABLE {}(classified TEXT, readname TEXT PRIMARY KEY, taxid INT)"""
        for name in tables:
            self.db.execute(drop_table.format(name))
            self.db.execute(create_table.format(name))

        self.reads = {}
        self.shared = defaultdict(dict)


    @staticmethod
    def parse_kaiju(kaiju_fn):
        with open(kaiju_fn) as f:
            for line in f:
                splitline = line.split()
                classification = splitline[0]
                readname = splitline[1]
                taxid = int(splitline[2])
                yield classification, readname, taxid 

    @staticmethod
    def parse_kraken(kraken_fn):
        with open(kraken_fn) as f:
            for line in f:
                splitline = line.split()
                classification = splitline[0]
                readname = splitline[1]
                taxid = int(splitline[2])
                yield classification, readname, taxid 

    @staticmethod
    def parse_clarks(clarks_fn):
        with open(clarks_fn) as f:
            f.readline() # Skip header line
            for line in f:
                splitline = line.split(",")
                classification = "C"
                if splitline[0].endswith(("/1", "/2")):
                    readname = splitline[0][:-2]
                else:
                    readname = splitline[0]
                try:
                    taxid = int(splitline[2])
                except ValueError:
                    taxid = 0
                    classification = "U"
                except IndexError:
                    print(line)
                yield classification, readname, taxid

    def fill_table(self, tablename, file_parser, lines_per_chunk=100000):
        """Inserts data into SQLite3 table.

        Implements some wild ideas on performance optimization from:
        http://stackoverflow.com/questions/1711631/improve-insert-per-second-performance-of-sqlite
        """
        logger.debug("Reading %s...", tablename)
        tic = time.time()
        insert_cmd = """INSERT INTO {table} VALUES (?, ?, ?)""".format(table=tablename)
        for chunk in grouper(lines_per_chunk, file_parser):
            self.db.executemany(insert_cmd, chunk)
        self.db.commit()
        toc = time.time()
        logger.debug("Reading %s completed in %2.2f seconds.", tablename, toc-tic)

    def compute_read_counts(self):
        logger.debug("Computing read counts for all data sources...")
        self.reads["kaiju"] = int(self.db.execute("SELECT Count(*) FROM kaiju").fetchone()[0])
        self.reads["kraken"] = int(self.db.execute("SELECT Count(*) FROM kraken").fetchone()[0])
        self.reads["clarks"] = int(self.db.execute("SELECT Count(*) FROM clarks").fetchone()[0])

        self.reads["classified_kaiju"] = int(self.db.execute('SELECT Count(classified) FROM kaiju WHERE classified = "C"').fetchone()[0])
        self.reads["classified_kraken"] = int(self.db.execute('SELECT Count(classified) FROM kraken WHERE classified = "C"').fetchone()[0])
        self.reads["classified_clarks"] = int(self.db.execute('SELECT Count(classified) FROM clarks WHERE classified = "C"').fetchone()[0])

        self.reads["unclassified_kaiju"] = int(self.db.execute('SELECT Count(classified) FROM kaiju WHERE classified = "U"').fetchone()[0])
        self.reads["unclassified_kraken"] = int(self.db.execute('SELECT Count(classified) FROM kraken WHERE classified = "U"').fetchone()[0])
        self.reads["unclassified_clarks"] = int(self.db.execute('SELECT Count(classified) FROM clarks WHERE classified = "U"').fetchone()[0])

        count_unique_cmd = """SELECT Count(readname) FROM 
                                (SELECT readname FROM kaiju
                                 UNION
                                 SELECT readname FROM kraken
                                 UNION 
                                 SELECT readname FROM clarks
                                )
                            """
        self.reads["total"] = int(self.db.execute(count_unique_cmd).fetchone()[0])

    def count_overlaps(self):
        logger.debug("Computing overlaps...")
        cmd_classified = """SELECT Count(*) FROM
                                (SELECT readname FROM {} WHERE classified = "C"
                                 INTERSECT
                                 SELECT readname FROM {} WHERE classified = "C")
                         """
        cmd_unclassified = """SELECT Count(*) FROM
                                 (SELECT readname FROM {} WHERE classified = "U"
                                  UNION
                                  SELECT readname FROM {} WHERE classified = "U")
                          """
        sources = ["kaiju", "kraken", "clarks"]
        Stuple = namedtuple("shared", "classified unclassified")
        for combination in sorted(chain(combinations(sources, 2))):
            classified = int(self.db.execute(cmd_classified.format(*combination)).fetchone()[0])
            unclassified = int(self.db.execute(cmd_unclassified.format(*combination)).fetchone()[0])
            self.shared[combination[0]][combination[1]] = Stuple(classified, unclassified)
            self.shared[combination[1]][combination[0]] = Stuple(classified, unclassified)

            if classified + unclassified != self.reads[combination[0]]:
                logger.warning("%s and %s are incongruent;  C:%s and U:%s should sum to %s",
                        combination[0], combination[1], classified, unclassified, self.reads[combination[0]])

    @staticmethod
    def validate_merge_order(merge_order):
        valid_names = {"kaiju", "kraken", "clarks"}
        for name in merge_order.split(","):
            if not name in valid_names:
                logger.error("Invalid name in merge order: %s", name)
                logger.error("Requested merge order: %s", merge_order)
                exit(2)
        logger.info("Merge order: %s", merge_order)

    def merge_tables(self, merge_order):
        logger.info("Merging data sources according to merge order: %s", merge_order)
        merge_order = merge_order.split(",")
        insert_first = """INSERT INTO merged SELECT * FROM {table}"""
        insert_or_replace = """INSERT OR REPLACE INTO merged 
                               SELECT * FROM {table} WHERE {table}.classified = "C"
                            """
        self.db.execute(insert_first.format(table=merge_order[0]))
        self.db.execute(insert_or_replace.format(table=merge_order[1]))
        self.db.execute(insert_or_replace.format(table=merge_order[2]))

        logger.debug("%s reads in merged table.", 
                self.db.execute("SELECT Count(readname) FROM merged;").fetchone()[0])

    def get_merged(self):
        for row in self.db.execute("SELECT * FROM merged ORDER BY readname"):
            yield map(str, row)


def main(dbfile=":memory:", output_fn="", kaiju="", kraken="", clarks="", merge_order=""):
    merger = Merger(dbfile)
    merger.validate_merge_order(merge_order)

    if kaiju:
        merger.fill_table("kaiju", merger.parse_kaiju(kaiju))
    if kraken:
        merger.fill_table("kraken", merger.parse_kraken(kraken))
    if clarks:
        merger.fill_table("clarks", merger.parse_clarks(clarks))

    merger.compute_read_counts()
    logger.info("Unique reads:")
    logger.info("  Kaiju:    %s", merger.reads["kaiju"])
    logger.info("  Kraken:   %s", merger.reads["kraken"])
    logger.info("  CLARK-S:  %s", merger.reads["clarks"])
    logger.info("  Combined: %s", merger.reads["total"])
    logger.info("Classified reads:")
    logger.info("  Kaiju:    %s", merger.reads["classified_kaiju"])
    logger.info("  Kraken:   %s", merger.reads["classified_kraken"])
    logger.info("  CLARK-S:  %s", merger.reads["classified_clarks"])
    logger.info("Unclassified reads:")
    logger.info("  Kaiju:    %s", merger.reads["unclassified_kaiju"])
    logger.info("  Kraken:   %s", merger.reads["unclassified_kraken"])
    logger.info("  CLARK-S:  %s", merger.reads["unclassified_clarks"])

    merger.merge_tables(merge_order)
    merger.count_overlaps()

    if output_fn:
        output = open(output_fn, 'w')
    else:
        output = stdout
    for row in merger.get_merged():
        print("\t".join(row), file=output)
    output.close()

    

if __name__ == "__main__":
    options, logger = parse_args()
    main(dbfile=options.dbfile,
            output_fn=options.output,
            kaiju=options.kaiju,
            kraken=options.kraken,
            clarks=options.clarks,
            merge_order=options.merge_order)

    
