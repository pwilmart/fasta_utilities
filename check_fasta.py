"""'check_fasta.py' written by PHil Wilmarth, OHSU, 2011-2020
# written by Phil Wilmarth, OHSU, 2011, 2020

The MIT License (MIT)

Copyright (c) 2020 Phillip A. Wilmarth and OHSU

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
# simple FASTA checking program - counts unusual amino acid characters

import os
import sys
import copy
import time
import fasta_lib

def fasta_checker(fasta_file, write):
    """Checks FASTA files for non-standard amino acid characters.
    """
    for obj in write:
        print("  database:", os.path.basename(fasta_file), file=obj)

    # initializations
    proteins = []
    p = fasta_lib.Protein()

    # counters
    prot_count = 0
    no_start_met = 0
    stop_count = 0
    stop_end = 0
    gap_count = 0
    B_count = 0
    J_count = 0
    O_count = 0
    U_count = 0
    X_count = 0
    Z_count = 0

    # read the sequences into a list
    f = fasta_lib.FastaReader(fasta_file)
    while f.readNextProtein(p, check_for_errs=True):
        prot_count += 1

        # test for odd amino acids, stop codons, gaps
        if not p.sequence.startswith('M'):
            no_start_met += 1
        if p.sequence.endswith('*'):
            stop_end += 1
        if '*' in p.sequence:
            stop_count += 1
        if '-' in p.sequence:
            gap_count += 1
        if 'B' in p.sequence:
            B_count += 1
        if 'J' in p.sequence:
            J_count += 1
        if 'O' in p.sequence:
            O_count += 1
        if 'U' in p.sequence:
            U_count += 1
        if 'X' in p.sequence:
            X_count += 1
        if 'Z' in p.sequence:
            Z_count += 1

        # save the protein in list
        proteins.append(copy.deepcopy(p))

    # check for duplicates and count
    duplicate_count = 0
    mw_dict = {}
    for i, p in enumerate(proteins):
        if mw_dict.get(str(p.molwtProtein()), False):
            j = mw_dict[str(p.molwtProtein())]
            if p.sequence == proteins[j].sequence:
                duplicate_count += 1
        else:
            mw_dict[str(p.molwtProtein())] = i

    # print out the report of oddball characters
    for obj in write:
        print("  total number of input sequences was:", prot_count, file=obj)
        print("  number of redundant sequences was:", duplicate_count, file=obj)
        print("    translations that do not start with Met:", no_start_met, file=obj)
        print("    translations that ended with a stop codon:", stop_end, file=obj)
        print("    translations that had premature stop codons:", stop_count, file=obj)
        print("    translations that contained gaps:", gap_count, file=obj)
        print("    translations that had B (ambiguous N/D):", B_count, file=obj)
        print("    translations that had J (ambiguous I/L):", J_count, file=obj)
        print("    translations that had O (pyrrolysine):", O_count, file=obj)
        print("    translations that had U (selenocysteine):", U_count, file=obj)
        print("    translations that had X (unknown amino acid):", X_count, file=obj)
        print("    translations that had Z (ambiguous Q/E):", Z_count, file=obj)

    return


# print program name and version
print('==========================================================')
print(' program check_fasta.py, v1.0.0, Phil Wilmarth, OHSU 2020 ')
print('==========================================================')

# browse to the database
database = r"C:\Xcalibur\database"
if not os.path.exists(database):
    database = os.getcwd()
file_ext_list = [('FASTA files', '*.fasta'), ('FASTA files', '*.fa'),
                 ('FASTA files', '*.gz')]
fasta_files = fasta_lib.get_files(database, file_ext_list, 'Select a FASTA database')
if not fasta_files:
    sys.exit()     # cancel button repsonse

# create a log file to mirror screen output
_folder = os.path.split(fasta_files[0])[0]
log_obj = open(os.path.join(_folder, 'fasta_utilities.log'), 'a')
write = [None, log_obj]
fasta_lib.time_stamp_logfile('\n>>> starting: check_fasta.py', log_obj)

# process the FASTA files
for fasta_file in fasta_files:
    try:
        fasta_checker(fasta_file, write)
    except FileNotFoundError:
        pass
    for obj in write:
        print(file=obj)

# finish up the log file
fasta_lib.time_stamp_logfile('>>> ending: check_fasta.py', log_obj)
log_obj.close()

# end
