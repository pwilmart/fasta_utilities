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
# simple TriTryp FASTA protein parsing/repair program
# truncates at stop codons, flags other odd things
# reformats description strings
# updated for Python 3 -PW 7/6/2017

import os
import copy
import fasta_lib_Py3 as fasta_lib

# print program name and version
print('============================================================')
print(' program TriTryp_fixer.py, v1.0.2, Phil Wilmarth, OHSU 2017 ')
print('============================================================')

# browse to the database
database = r"C:\Xcalibur\database"
if not os.path.exists(database):
    database = os.getcwd()
fasta_file = fasta_lib.get_file(database, [('FASTA files', '*.fasta')], 'Select a TriTryp FASTA database')
if fasta_file == '': sys.exit()     # cancel button repsonse

# build new database name
new_fasta_file = os.path.basename(fasta_file)
new_fasta_file = new_fasta_file.replace('.fasta', '_fixed.fasta')
new_fasta_file = os.path.join(os.path.dirname(fasta_file), new_fasta_file)

# initializations
proteins = []
p = fasta_lib.Protein()
pcount = 0
stop_count = 0
gap_count = 0
no_met = 0

# read the sequences into a list
f = fasta_lib.FastaReader(fasta_file)
while f.readNextProtein(p, check_for_errs=True):
    pcount += 1
    
    # parse the description string into a dictionary
    items = [x.strip() for x in p.description.split('|') if x]
    header_dict = {x.split('=')[0]: x.split('=')[1] for x in items}
    new_desc = []
    new_desc = [header_dict['gene_product'],
                '(' + header_dict['location'] + ')',
                '[' + header_dict['organism'] + ']',
                '(' + header_dict['protein_length'] + 'aa)']
    p.new_desc = ' '.join(new_desc)
    
    # test for odd amino acids, stop codons, gaps
    if not p.sequence.startswith('M'):
        no_met += 1
        p.new_desc = p.new_desc + ' (No starting Met)'
    if '*' in p.sequence:
        stop_count += 1
        cut = p.sequence.index('*')
        string = ' (Premature stop %s/%s)' % (cut, len(p.sequence))
        p.new_desc = p.new_desc + string
        p.sequence = p.sequence[:cut]
    if '-' in p.sequence:
        gap_count += 1
        p.new_desc = p.new_desc + ' (Contains gaps)'
    
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

# open the new protein fasta file and write out the proteins
file_obj = open(new_fasta_file, 'w')
for p in proteins:
    p.printProtein(file_obj)
file_obj.close()

# print out the report of oddball characters
print("   TriTryp database:", os.path.basename(fasta_file))
print("   translations that do not start with Met:", no_met)
print("   translations that have premature stop codons:", stop_count)
print("   translations that contain gaps:", gap_count)
print("   total number of input sequences was:", pcount)
print("   total number of sequences written was:", len(proteins))
print("   number of redundant sequences was:", duplicates)

# end
