"""'reverse_fasta.py' Written by Phil Wilmarth, OHSU.

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
# updated for Python 3 -PW 7/7/2017

import os
import sys
import fasta_lib_Py3 as fasta_lib

# flag to parse (clean) accessions and descriptions.
# NOTE: cleaning can cause some loss of information.  Apply at later stages
# of database processing such as here.
CLEAN_ACCESSIONS = False

# flag to keep UniProt Identifier (True) or more-stable ACC (False)
KEEP_UNIPROT_ID = False

# flag set if RefSeq entries are being extracted from NCBI nr
REF_SEQ_ONLY = True

# flags to make different output databases
MAKE_SEPARATE_FORWARD = True
MAKE_SEPARATE_REVERSED = False
MAKE_SEPARATE_BOTH = True


def fasta_reverse(fasta_file):
    """Adds contaminants and reverses entries in a FASTA protein database.

    Called with FASTA filename.  Reversed DB written to same location.
    Options for separate or concatenated output files.
    """
    decoy_string = 'REV_'   # the string to denote decoy sequences
    
    # open the "forward" and "reversed" output files
    if fasta_file.lower().endswith('.gz'):
        _file = os.path.splitext(fasta_file[:-3])[0]
    else:
        _file = os.path.splitext(fasta_file)[0]
    for_name = _file + '_for.fasta'
    for_file_obj = open(for_name, 'w')
    rev_name = _file + '_rev.fasta'
    rev_file_obj = open(rev_name, 'w')
    
    # create a log file to mirror screen output
    _folder = os.path.split(fasta_file)[0]
    log_obj = open(os.path.join(_folder, 'fasta_utilities.log'), 'a')
    write = [None, log_obj]
    fasta_lib.time_stamp_logfile('\n>>> starting: reverse_fasta.py', log_obj)
    
    # create instances protein object and initialize counter
    prot = fasta_lib.Protein()
    pcount = 0

    # try to find the contaminants database file
    try:
        if os.path.exists('all_contams_fixed.fasta'):
            print('contams file found:', _file)
            _file = 'all_contams_fixed.fasta'
        else:
            path = os.path.split(fasta_file)[0]
            _file = os.path.join(path, 'all_contams_fixed.fasta')
            print('trying:', _file)

        # create reader and fetch contaminants           
        f = fasta_lib.FastaReader(_file)
        while f.readNextProtein(prot, check_for_errs=True):
            pcount += 1
            if CLEAN_ACCESSIONS:
                prot.parseCONT()
            prot.printProtein(for_file_obj)
            rev = prot.reverseProtein(decoy_string)
            rev.printProtein(rev_file_obj)
        for obj in write:
            print('...there were %s contaminant entries in %s' %
                  ("{0:,d}".format(pcount), _file), file=obj)           
    except:
        for obj in write:
            print('   WARNING: "all_contams_fixed.fasta" not found!', file=obj)
    
    # read proteins until EOF and write proteins to "forward" and "reversed" files
    f = fasta_lib.FastaReader(fasta_file)
    
    # error checking slows program execution, turn on if needed.
    # Reading and writing sequences always removes spaces and blank lines.
    while f.readNextProtein(prot, check_for_errs=False):
        pcount += 1
        if CLEAN_ACCESSIONS:
            if prot.accession.startswith('gi|'):
                prot.parseNCBI(REF_SEQ_ONLY)
            elif prot.accession.startswith('sp|') or prot.accession.startswith('tr|'):
                prot.parseUniProt(KEEP_UNIPROT_ID)
        
        prot.printProtein(for_file_obj)    # write to "forward" file
        rev = prot.reverseProtein(decoy_string)
        rev.printProtein(rev_file_obj)   # write to "reversed" file
    
    # make concatenated output file if desired and print summary stats
    if MAKE_SEPARATE_BOTH:
        _file = fasta_file.replace('.gz', '')
        both_name = _file.replace('.fasta', '_both.fasta')
        both_file_obj = open(both_name, 'w')
        for_file_obj.close()
        for_file_obj = open(for_name, 'r')
        rev_file_obj.close()
        rev_file_obj = open(rev_name, 'r')
        while True:
            line = for_file_obj.readline()
            if not line: break
            both_file_obj.write(str(line))
        while True:
            line = rev_file_obj.readline()
            if not line: break
            both_file_obj.write(str(line))
        both_file_obj.close()
        for obj in write:
            print('...%s total proteins written to %s' %
                  ("{0:,d}".format(2*pcount), os.path.split(both_name)[1]), file=obj)
    
    if MAKE_SEPARATE_FORWARD:
        for obj in write:
            print('...%s proteins written to %s' %
                  ("{0:,d}".format(pcount), os.path.split(for_name)[1]), file=obj)
    if MAKE_SEPARATE_REVERSED:
        for obj in write:
            print('...%s proteins reversed and written to %s' %
                  ("{0:,d}".format(pcount), os.path.split(rev_name)[1]), file=obj)
    
    # close files and delete unwanted files
    for_file_obj.close()
    rev_file_obj.close()
    fasta_lib.time_stamp_logfile('>>> ending: reverse_fasta.py', log_obj)
    log_obj.close()
    if not MAKE_SEPARATE_FORWARD:
        os.remove(for_name)
    if not MAKE_SEPARATE_REVERSED:
        os.remove(rev_name)
    return


# setup stuff: check for command line args, etc.
if __name__ == '__main__':
    # check if database name passed on command line
    if len(sys.argv) > 1:
        fasta_file = sys.argv[1:]
    
    # if not, browse to database file
    else:
        if len(sys.argv) > 1:
            print('...WARNING: %s not found...' % (sys.argv[1],))
        database = r'C:\Xcalibur\database'
        if not os.path.exists(database):
            database = os.getcwd()
        fasta_files = fasta_lib.get_files(database,
                                          [('FASTA files', '*.fasta'), ('Zipped files', '*.gz'), ('All files', '*.*')],
                                          'Select a FASTA database')
        if not fasta_files: sys.exit() # cancel button response

    # print version info, etc. (here because of the loop)    
    print('=================================================================')
    print(' reverse_fasta.py, v 1.1.0, written by Phil Wilmarth, OHSU, 2017 ')
    print('=================================================================')
    
    # call reverse function
    os.chdir('.')   # set location to where script lives - the contaminats should be there
    for fasta_file in fasta_files:
        try:
            fasta_reverse(fasta_file)
        except FileNotFoundError:   # FastaReader class raises exception if file not found
            pass

# end
