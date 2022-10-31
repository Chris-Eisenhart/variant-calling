#!/usr/bin/env python3.9

"""
Toolkit entry point for variant calling.
"""
import argparse
import os
import sys

from variant_calling import call_variants_on_sam_file


def parseArgs(args): 
    """
    Parse the command line arguments into useful python objects.  '--' variables are optional
    set_defaults only applies when the argument is not provided (it won't override)
    """
    parser = argparse.ArgumentParser(description = __doc__)
    parser.add_argument("sam_file",
                        help = " The input sam file",
                        action = "store")
    parser.add_argument("out_variant_file",
                        help = " The output variant file",
                        action = "store")
    parser.add_argument("out_coverage_file",
                        help = " The output coverage file",
                        action = "store")
    parser.add_argument("--verbose",
                        help = " The verbosity level for stdout messages (default INFO)",
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                        action = "store")
    parser.set_defaults(verbose = "INFO")
    options = parser.parse_args()
    return options


def main(args):
    options = parseArgs(args)
    call_variants_on_sam_file(options.sam_file, options.out_variant_file, options.out_coverage_file)


if __name__ == "__main__" :
    sys.exit(main(sys.argv))
