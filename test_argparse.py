"""test_argparse.py
written by Phil Wilmarth, 2017.
"""
import os
import sys
import argparse

MAKE_FORWARD = False
MAKE_REVERSE = False
MAKE_BOTH = True

parser = argparse.ArgumentParser(description='Makes databases with contaminants and decoys.',
                                 prefix_chars='-+')
parser.add_argument('+f', '++forward', dest='forward',
                    help='makes forward sequences with contaminants',
                    action='store_true', default=MAKE_FORWARD)
parser.add_argument('-f', '--forward', dest='forward',
                    help='does not makes forward sequences with contaminants',
                    action='store_false', default=MAKE_FORWARD)
parser.add_argument('+r', '++reverse', dest='reverse',
                    help='makes reversed sequences with contaminants',
                    action='store_true', default=MAKE_REVERSE)
parser.add_argument('-r', '--reverse', dest='reverse',
                    help='does not makes reversed sequences with contaminants',
                    action='store_false', default=MAKE_REVERSE)
parser.add_argument('+b', '++both', dest='both',
                    help='makes forward and reversed sequences with contaminants',
                    action='store_true', default=MAKE_BOTH)
parser.add_argument('-b', '--both', dest='both',
                    help='does not makes forward and reversed sequences with contaminants',
                    action='store_false', default=MAKE_BOTH)
parser.add_argument('-v', '--version', action='version', version='%(prog)s version 1.1.1')
parser.add_argument('file', help='list of FASTA files to process', nargs='*')

args = parser.parse_args()

print('forward:', args.forward)
print('reverse:', args.reverse)
print('both:', args.both)
print('file(s):', args.file)

print('args:', args)
print(sys.argv)
