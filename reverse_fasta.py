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
# improved command line agrument processing -PW 7/30/2017

import os
import sys
import argparse
import fasta_lib

# flags to make different output databases
MAKE_FORWARD = True
MAKE_REVERSE = False
MAKE_BOTH = True


def main(fasta_file, forward=False, reverse=False, both=True, log_obj=None, contam_path=""):
    """Adds contaminants and reverses entries for a FASTA protein database.

    Call with single fasta file name.
    If "forward", make sequences plus contaminants,
    if "reverse", make reversed sequences with reversed contaminants,
    if "both", make concatenated target/decoy with contaminants.
    "contam_path" is optional fullpath name of a contaminants database to use instead of default
    """
    decoy_string = 'REV_'   # the string to denote decoy sequences
    ######################################
    # Change default contaminants file name here:
    CONTAMS = 'Thermo_contams_fixed.fasta'
    # or pass in a "contams_path"
    ######################################
    
    # open the "forward" and "reversed" output files
    if fasta_file.lower().endswith('.gz'):
        _file = os.path.splitext(fasta_file[:-3])[0]
    else:
        _file = os.path.splitext(fasta_file)[0]
    for_name = _file + '_for.fasta'
    for_file_obj = open(for_name, 'w')
    rev_name = _file + '_rev.fasta'
    rev_file_obj = open(rev_name, 'w')

    # create the name for the concatenated file (if later needed)
    both_name = _file + '_both.fasta'
    
    # create a log file to mirror screen output
    _folder = os.path.split(fasta_file)[0]
    if not log_obj:
        log_obj = open(os.path.join(_folder, 'fasta_utilities.log'), 'a')
    write = [None, log_obj]
    fasta_lib.time_stamp_logfile('\n>>> starting: reverse_fasta.py', log_obj)
    
    # create instances protein object and initialize counter
    prot = fasta_lib.Protein()
    pcount = 0

    # try to find the contaminants database file
    # If no contam file path provided, search for it in current directory
    _file = None
    if not contam_path:
        if os.path.exists(CONTAMS):
            _file = CONTAMS
        else:
            path = os.path.split(fasta_file)[0]
            if os.path.exists(os.path.join(path, CONTAMS)):
                _file = os.path.join(path, CONTAMS)
    elif os.path.exists(contam_path) and os.path.isfile(contam_path):
        _file = contam_path
    elif os.path.isdir(contam_path) and os.path.exists(os.path.join(contam_path, CONTAMS)):
        _file = os.path.join(contam_path, CONTAMS)
        
    # create reader and add contaminants (if contams file was found)
    if _file:
        f = fasta_lib.FastaReader(_file)
        while f.readNextProtein(prot, check_for_errs=True):
            pcount += 1
            prot.printProtein(for_file_obj)
            rev = prot.reverseProtein(decoy_string)
            rev.printProtein(rev_file_obj)
        for obj in write:
            print('...there were %s contaminant entries in %s' %
                  ("{0:,d}".format(pcount), os.path.split(_file)[1]), file=obj)
    else:        
        for obj in write:
            print('...WARNING: contaminants were not added', file=obj)
        
    # read proteins until EOF and write proteins to "forward" and "reversed" files
    f = fasta_lib.FastaReader(fasta_file)
    
    # error checking slows program execution, turn on if needed.
    # Reading and writing sequences always removes spaces and blank lines.
    while f.readNextProtein(prot, check_for_errs=False):
        pcount += 1
        prot.printProtein(for_file_obj)    # write to "forward" file
        rev = prot.reverseProtein(decoy_string)
        rev.printProtein(rev_file_obj)   # write to "reversed" file
    for_file_obj.close()
    rev_file_obj.close()
    
    # make concatenated output file if desired and print summary stats
    if both:
        both_file_obj = open(both_name, 'w')
        for_file_obj = open(for_name, 'r')
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
    
    if forward:
        for obj in write:
            print('...%s proteins written to %s' %
                  ("{0:,d}".format(pcount), os.path.split(for_name)[1]), file=obj)
    if reverse:
        for obj in write:
            print('...%s proteins reversed and written to %s' %
                  ("{0:,d}".format(pcount), os.path.split(rev_name)[1]), file=obj)
    
    # close files and delete unwanted files
    for_file_obj.close()
    rev_file_obj.close()
    fasta_lib.time_stamp_logfile('>>> ending: reverse_fasta.py', log_obj)
    log_obj.close()
    if not forward:
        os.remove(for_name)
    if not reverse:
        os.remove(rev_name)
    return


# setup stuff: check for command line args, etc.
if __name__ == '__main__':
    # set up command line arguments
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
    parser.add_argument('-v', '--version', action='version', version='%(prog)s version 1.1.2')
    parser.add_argument('files', help='list of FASTA files to process', nargs='*')

    args = parser.parse_args()
        
    if len(sys.argv) > 1:   # options set from command line
        forward = args.forward
        reverse = args.reverse
        both = args.both
        fasta_files = args.files
    else:   # options set to hardcoded defaults if interactive mode or no passed commands
        forward = MAKE_FORWARD
        reverse = MAKE_REVERSE
        both = MAKE_BOTH
        fasta_files = []
    
    # if no FASTA files, browse to database file(s)
    if not fasta_files:
        database = r'C:\Xcalibur\database'
        if not os.path.exists(database):
            database = os.getcwd()
        fasta_files = fasta_lib.get_files(database,
                                          [('FASTA files', '*.fasta'), ('Zipped files', '*.gz'), ('All files', '*.*')],
                                          'Select a FASTA database')
        if not fasta_files: sys.exit() # cancel button response

    # print version info, etc. (here because of the loop)    
    print('=================================================================')
    print(' reverse_fasta.py, v 1.1.2, written by Phil Wilmarth, OHSU, 2017 ')
    print('=================================================================')
    
    # call reverse function
    os.chdir('.')   # set location to where script lives - the contaminants should be there
    for fasta_file in fasta_files:
        try:
            main(fasta_file, forward, reverse, both)
        except IOError:   # FastaReader class raises exception if file not found
            print('...WARNING: %s not found' % fasta_file)
            pass

# end
