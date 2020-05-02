"""'TryTrip_fixer.py' written by PHil Wilmarth, OHSU, 2011-2017
# written by Phil Wilmarth, OHSU, 2011

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
# simple FASTA checking program - counts unusual things

import os
import copy
import fasta_lib

# print program name and version
print('============================================================')
print(' program FASTA_checker.py, v1.0.0, Phil Wilmarth, OHSU 2020 ')
print('============================================================')

# browse to the database
database = r"C:\Xcalibur\database"
if not os.path.exists(database):
    database = os.getcwd()
fasta_file = fasta_lib.get_file(database, [('FASTA files', '*.fasta')], 'Select a FASTA database')
if fasta_file == '': sys.exit()     # cancel button repsonse

# initializations
proteins = []
p = fasta_lib.Protein()
pcount = 0
stop_count = 0
stop_end = 0
gap_count = 0
no_met = 0

# read the sequences into a list
f = fasta_lib.FastaReader(fasta_file)
while f.readNextProtein(p, check_for_errs=True):
    pcount += 1

    # test for odd amino acids, stop codons, gaps
    if not p.sequence.startswith('M'):
        no_met += 1
    if p.sequence.endswith('*'):
        stop_end += 1
    if '*' in p.sequence:
        stop_count += 1
    if '-' in p.sequence:
        gap_count += 1

    # save the protein in list
    proteins.append(copy.deepcopy(p))

# check for duplicates and count
duplicates = 0
mw_dict = {}
for i, p in enumerate(proteins):
    if mw_dict.get(str(p.molwtProtein()), False):
        j = mw_dict[str(p.molwtProtein())]
        if p.sequence == proteins[j].sequence:
            duplicates += 1
    else:
        mw_dict[str(p.molwtProtein())] = i

# print out the report of oddball characters
print("   database:", os.path.basename(fasta_file))
print("   translations that do not start with Met:", no_met)
print("   translations that have premature stop codons:", stop_count)
print("   translations that ended with a stop codon:", stop_end)
print("   translations that contain gaps:", gap_count)
print("   total number of input sequences was:", pcount)
print("   number of redundant sequences was:", duplicates)

# end
