#!/usr/bin/env python3
__author__ = "Fredrik Boulund"
__date__ = "2016"
__doc__ = """Merge taxonomic read assignments produced by Kaiju, Kraken, and CLARK-S."""
__version__ = "1.0"

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

    desc = """Merge read classifications from Kaiju, Kraken, and CLARK-S. 
              Copyright (c) {author} {date}. 
              Version {version}.""".format(author=__author__, 
                                           date=__date__,
                                           version=__version__)
    epilog = """Redefine merge order to decide classification priority; 
    classifications are overwritten according to merge order.
    Using an sqlite3 database on disk consumes <150MB RAM. This is recommended.
    The database consumes about the sum of the sizes of the input files."""

    parser = argparse.ArgumentParser(description=desc, epilog=epilog)

    parser.add_argument("-k", "--kaiju", dest="kaiju", metavar="KAIJU",
            help="Kaiju output file (tab).")
    parser.add_argument("-K", "--kraken", dest="kraken", metavar="KRAKEN",
            help="Kraken output file (tab).")
    parser.add_argument("-c", "--clarks", dest="clarks", metavar="CLARKS",
            help="CLARK-S output file (csv).")
    parser.add_argument("-m", "--merge-order", dest="merge_order", metavar="ORDER",
            default="clarks,kraken,kaiju",
            help="Define merge order, e.g. kraken,kaiju,clarks [%(default)s].")
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
    parser.add_argument("--version", action="store_true",
            default=False,
            help="Print version.")

    if len(argv) < 2:
        parser.print_help()
        exit(1)

    options = parser.parse_args()

    if options.version:
        print(__doc__+" Version {}".format(__version__))
        exit()

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
        
        # Optimizations for fast inserts (see self.insert_first)
        self.db.execute("PRAGMA journal_mode = MEMORY")
        self.db.execute("PRAGMA synchronous = OFF")
        
        drop_table = """DROP TABLE IF EXISTS merged"""
        create_table = """CREATE TABLE merged (
                             classified TEXT, 
                             readname TEXT PRIMARY KEY, 
                             taxid INT,
                             source TEXT
                             )"""
        self.db.execute(drop_table)
        self.db.execute(create_table)

        self.reads = {}
        self.shared = defaultdict(dict)

    def validate_merge_order(self, merge_order):
        valid_names = {"kaiju": ("kaiju", self.parse_kaiju),
                       "kraken": ("kraken", self.parse_kraken),
                       "clarks": ("clarks", self.parse_clarks)}
        valid_merge_order = []
        for name in merge_order.split(","):
            try:
                valid_merge_order.append(valid_names[name])
            except KeyError:
                logger.error("Invalid name in merge order: %s", name)
                logger.error("Requested merge order: %s", merge_order)
                exit(2)
        logger.info("Merge order: %s", merge_order)
        return valid_merge_order

    def parse_kaiju(self, kaiju_fn, classified_only=False):
        with open(kaiju_fn) as f:
            classifications = 0
            for num, line in enumerate(f, start=1):
                splitline = line.split()
                classification = splitline[0]
                if classification == "C":
                    classifications += 1
                readname = splitline[1]
                taxid = int(splitline[2])
                if classified_only and classification == "C":
                    yield classification, readname, taxid, "kaiju"
                elif classified_only:
                    pass
                else:
                    yield classification, readname, taxid, "kaiju"
        self.reads["kaiju"] = num
        self.reads["classified_kaiju"] = classifications
        self.reads["unclassified_kaiju"] = num - classifications

    def parse_kraken(self, kraken_fn, classified_only=False):
        with open(kraken_fn) as f:
            classifications = 0
            for num, line in enumerate(f, start=1):
                splitline = line.split()
                classification = splitline[0]
                if classification == "C":
                    classifications += 1
                readname = splitline[1]
                taxid = int(splitline[2])
                if classified_only and classification == "C":
                    yield classification, readname, taxid, "kraken"
                elif classified_only:
                    pass
                else:
                    yield classification, readname, taxid, "kraken"
        self.reads["kraken"] = num
        self.reads["classified_kraken"] = classifications
        self.reads["unclassified_kraken"] = num - classifications

    def parse_clarks(self, clarks_fn, classified_only=False):
        with open(clarks_fn) as f:
            f.readline() # Skip header line
            classifications = 0
            for num, line in enumerate(f, start=1):
                splitline = line.split(",")
                if splitline[0].endswith(("/1", "/2")):
                    readname = splitline[0][:-2]
                else:
                    readname = splitline[0]
                try:
                    taxid = int(splitline[2])
                    classification = "C"
                    classifications += 1
                except ValueError:
                    taxid = 0
                    classification = "U"
                except IndexError:
                    logger.error("Could not parse CLARK-S line:")
                    logger.error(line)
                if classified_only and classification == "C":
                    yield classification, readname, taxid, "clarks"
                elif classified_only:
                    pass
                else:
                    yield classification, readname, taxid, "clarks"
        self.reads["clarks"] = num
        self.reads["classified_clarks"] = classifications
        self.reads["unclassified_clarks"] = num - classifications

    def insert_first(self, source_name, file_parser, lines_per_chunk=100000):
        """Inserts data into SQLite3 table.

        Implements some wild ideas on performance optimization from:
        http://stackoverflow.com/questions/1711631/improve-insert-per-second-performance-of-sqlite
        """
        logger.debug("Reading %s...", source_name)
        tic = time.time()
        insert_cmd = """INSERT INTO merged VALUES (?, ?, ?, ?)"""
        for chunk in grouper(lines_per_chunk, file_parser):
            self.db.executemany(insert_cmd, chunk)
        self.db.commit()
        toc = time.time()
        logger.debug("Reading %s completed in %2.2f seconds.", source_name, toc-tic)

    def count_sources(self):
        logger.debug("Summarizing merged classifications...")
        summary_cmd = """SELECT source, classified, Count(classified) FROM merged 
                                GROUP BY source, classified
                         """
        for row in self.db.execute(summary_cmd).fetchall():
            yield row
        logger.debug("Summarizing merged classifications completed.")

    def insert_replace(self, source_name, file_parser, lines_per_chunk=100000):
        logger.debug("Reading %s...", source_name)
        tic = time.time()
        insert_or_replace = """INSERT OR REPLACE INTO merged VALUES (?, ?, ?, ?)"""
        for chunk in grouper(lines_per_chunk, file_parser):
            self.db.executemany(insert_or_replace, chunk)
        self.db.commit()
        toc = time.time()
        logger.debug("Reading %s completed in %2.2f seconds.", source_name, toc-tic)

    def get_merged(self):
        for row in self.db.execute("SELECT classified, readname, taxid FROM merged ORDER BY readname"):
            yield map(str, row)


def main(dbfile=":memory:", output_fn="", kaiju="", kraken="", clarks="", merge_order=""):
    merger = Merger(dbfile)
    source_files = {"kaiju": kaiju,
                    "kraken": kraken,
                    "clarks": clarks}
    valid_merge_order = merger.validate_merge_order(merge_order)

    merger.insert_first(source_name=valid_merge_order[0][0],
                        file_parser=valid_merge_order[0][1](source_files[valid_merge_order[0][0]]))
    for source_name, parser in valid_merge_order[1:]:
        merger.insert_replace(source_name=source_name, 
                              file_parser=parser(source_files[source_name], 
                                                 classified_only=True))

    logger.info(" Source summary ".center(50, "="))
    logger.info("Number of reads:")
    for source_name, _ in valid_merge_order:
        logger.info("  %7s:  %s", source_name, merger.reads[source_name])

    logger.info("Classified reads:")
    for source_name, _ in valid_merge_order:
        logger.info("  %7s:    %9i (%2.2f%%)", 
                source_name, 
                merger.reads["classified_"+source_name], 
                100 * merger.reads["classified_"+source_name] / merger.reads[source_name])

    logger.info("Unclassified reads:")
    for source_name, _ in valid_merge_order:
        logger.info("  %7s:    %9i (%2.2f%%)", 
                source_name, 
                merger.reads["unclassified_"+source_name],
                100 * (merger.reads["unclassified_"+source_name]) / merger.reads[source_name])

    logger.info(" Merged classifications summary ".center(50, "="))
    counted_sources = list(merger.count_sources())
    classified_total = sum(c[2] for c in counted_sources if c[1] == "C")
    unclassified_total = sum(c[2] for c in counted_sources if c[1] == "U")
    if classified_total:
        logger.info(" Classified reads:")
    for source, count in ((row[0], row[2]) for row in counted_sources if row[1] == "C"):
        logger.info("  %7s: %10i", source, count)
    if unclassified_total:
        logger.info(" Unclassified reads:")
    for source, count in ((row[0], row[2]) for row in counted_sources if row[1] == "U"):
        logger.info("  %7s: %10i", source, count)
    logger.info(" Total classified:   %10i (%2.2f%%)", classified_total, 
            100 * classified_total/(classified_total+unclassified_total))
    logger.info(" Total unclassified: %10i (%2.2f%%)", unclassified_total, 
            100 * unclassified_total/(classified_total+unclassified_total))

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

    
