"""'count_fasta.py' Written by Phil Wilmarth, OHSU.

The MIT License (MIT)

Copyright (c) 2017 Phillip A. Wilmarth and OHSU

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

Direct questions to:
Technology & Research Collaborations, Oregon Health & Science University,
Ph: 503-494-8200, FAX: 503-494-4729, Email: techmgmt@ohsu.edu.
"""
# updated for Python 3 compatibility -PW 7/4/2017

import os
import sys
import time
import re
import fasta_lib


def fasta_counter(fasta_file):
    """Counts entries in a FASTA protein database.
        Call with FASTA filename.
        Checks for duplicate accessions and (optional) valid characters.
    """
    # create a log file to mirror screen output
    _folder = os.path.split(fasta_file)[0]
    log_obj = open(os.path.join(_folder, 'fasta_utilities.log'), 'a')
    write = [None, log_obj]
    fasta_lib.time_stamp_logfile('\n>>> starting: count_fasta.py', log_obj)

    # create instances of reader object and protein object, initialize counters
    f = fasta_lib.FastaReader(fasta_file)
    p = fasta_lib.Protein()
    prot = 0
    head = 0
    conflict = {}

    # construct summary file name based on FASTA name
    if fasta_file.endswith('.fasta.gz'):
        summary_file = re.sub(r'.fasta.gz$', r'.txt', fasta_file)
    else:
        summary_file = re.sub(r'.fasta$', r'.txt', fasta_file)
    if summary_file == fasta_file:
        summary_file = fasta_file + '.txt'

    # open summary file and write header
    summary_obj = open(summary_file, mode='wt')
    summary_obj.write('Accession\tLength\tMW\n')

    # read proteins until EOF; NOTE: checking for errors slows program by factor of 3-4
    while f.readNextProtein(p, check_for_errs=False):

        # count protein sequences
        prot += 1
        if (prot % 500000) == 0:
            print('......(%s proteins read...)' % ("{0:,d}".format(prot),))

##        # check for duplicate accession
##        dup = conflict.get(p.accession, False)
##        if dup:
##            for obj in write:
##                print('\n...WARNING: %s is already in FASTA database!\n' % (p.accession,), file=obj)
##                if p.molwtProtein(show_errs=False) == conflict[p.accession]:
##                    print('......possible duplicated sequence...', file=obj)
##        else:
##            conflict[p.accession] = p.molwtProtein(show_errs=False)

        # count number of header elements
        control_A = p.description.count(chr(1))
        head = head + control_A + 1

        # add info to summary_file
        print('\t'.join([p.accession, str(p.seqlenProtein()), str(round(p.molwtProtein(), 1))]), file=summary_obj)

    # print results and return
    for obj in write:
        print('...there are %s proteins in %s' % ("{0:,d}".format(prot), os.path.split(fasta_file)[1]), file=obj)
        if head > prot:
            print('...there were %s header lines' % ("{0:,d}".format(head),), file=obj)

    fasta_lib.time_stamp_logfile('>>> ending: count_fasta.py', log_obj)
    log_obj.close()
    summary_obj.close()
    return


# setup stuff: check for command line args, etc.
if __name__ == '__main__':
    # check if database name(s) passed on command line
    if len(sys.argv) > 1:
        fasta_files = sys.argv[1:]

    # if not, browse to database file
    else:
        database = r'C:\Xcalibur\database'  # set a default to speed up browsing
        if not os.path.exists(database):
            database = os.getcwd()
        fasta_files = fasta_lib.get_files(database,
                                          [('FASTA files', '*.fasta'), ('Zipped FASTA files', '*.gz'), ('All files', '*.*')],
                                          'Select a FASTA database')
        if not fasta_files: sys.exit()     # cancel button repsonse

    # print version info, etc. (here because of the loop)
    print('===============================================================')
    print(' count_fasta.py, v 1.1.0, written by Phil Wilmarth, OHSU, 2017 ')
    print('===============================================================')

    # call counter function for each fasta file
    print('start', time.ctime())
    for fasta_file in fasta_files:
        try:
            fasta_counter(fasta_file)
        except FileNotFoundError:     # FastaReader class raises exception if file not found
            pass
    print('end', time.ctime())

# end
