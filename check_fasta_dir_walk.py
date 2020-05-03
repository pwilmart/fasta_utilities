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

class Species():
    """generic data container."""
    def __init__(self):
        """Bare bones constructor."""
        self.fasta_file = None
        self.common_name = None
        self.latin_name = None
        self.taxon = None
        self.release = None
        self.db_name = None
        self.prot_count = None
        self.duplicate_count = None
        self.no_start_met_count = None
        self.stop_end_count = None
        self.stop_count = None
        self.gap_count = None
        self.B_count = None
        self.J_count = None
        self.O_count = None
        self.U_count = None
        self.X_count = None
        self.Z_count = None

    def make_header(self):
        headers = ['Release', 'FASTA database', 'Common name', 'Latin name', 'Taxonomy',
                   'Sequence count', 'Duplicate sequences', 'No start Met',
                   'Sequences with stop at end', 'Sequences with internal stop',
                   'Sequences with gap', 'Sequences with B', 'Sequences with J',
                   'Sequences with O', 'Sequences with U', 'Sequences with X',
                   'Sequences with Z']
        return '\t'.join(headers)

    def make_row(self):
        cells = [self.release, self.db_name, self.common_name, self.latin_name, self.taxon,
                 self.prot_count, self.duplicate_count, self.no_start_met_count,
                 self.stop_end_count, self.stop_count, self.gap_count, self.B_count,
                 self.J_count, self.O_count, self.U_count, self.X_count, self.Z_count]
        return '\t'.join([str(x) for x in cells])

def fasta_checker(species, write):
    """Checks FASTA files for non-standard amino acid characters.
    """
    for obj in write:
        print("  database:", os.path.basename(species.fasta_file), file=obj)

    # initializations
    proteins = []
    p = fasta_lib.Protein()

    # counters
    species.prot_count = 0
    species.no_start_met_count = 0
    species.stop_count = 0
    species.stop_end_count = 0
    species.gap_count = 0
    species.B_count = 0
    species.J_count = 0
    species.O_count = 0
    species.U_count = 0
    species.X_count = 0
    species.Z_count = 0

    # read the sequences into a list
    f = fasta_lib.FastaReader(species.fasta_file)
    while f.readNextProtein(p, check_for_errs=True):
        species.prot_count += 1

        # test for odd amino acids, stop codons, gaps
        if not p.sequence.startswith('M'):
            species.no_start_met_count += 1
        if p.sequence.endswith('*'):
            species.stop_end_count += 1
        if '*' in p.sequence:
            species.stop_count += 1
        if '-' in p.sequence:
            species.gap_count += 1
        if 'B' in p.sequence:
            species.B_count += 1
        if 'J' in p.sequence:
            species.J_count += 1
        if 'O' in p.sequence:
            species.O_count += 1
        if 'U' in p.sequence:
            species.U_count += 1
        if 'X' in p.sequence:
            species.X_count += 1
        if 'Z' in p.sequence:
            species.Z_count += 1

        # save the protein in list
        proteins.append(copy.deepcopy(p))

    # check for duplicates and count
    species.duplicate_count = 0
    mw_dict = {}
    for i, p in enumerate(proteins):
        if mw_dict.get(str(p.molwtProtein()), False):
            j = mw_dict[str(p.molwtProtein())]
            if p.sequence == proteins[j].sequence:
                species.duplicate_count += 1
        else:
            mw_dict[str(p.molwtProtein())] = i

    # print out the report of oddball characters
    for obj in write:
        print("  total number of input sequences was:", species.prot_count, file=obj)
        print("  number of redundant sequences was:", species.duplicate_count, file=obj)
        print("    translations that do not start with Met:", species.no_start_met_count, file=obj)
        print("    translations that ended with a stop codon:", species.stop_end_count, file=obj)
        print("    translations that had premature stop codons:", species.stop_count, file=obj)
        print("    translations that contained gaps:", species.gap_count, file=obj)
        print("    translations that had B (ambiguous N/D):", species.B_count, file=obj)
        print("    translations that had J (ambiguous I/L):", species.J_count, file=obj)
        print("    translations that had O (pyrrolysine):", species.O_count, file=obj)
        print("    translations that had U (selenocysteine):", species.U_count, file=obj)
        print("    translations that had X (unknown amino acid):", species.X_count, file=obj)
        print("    translations that had Z (ambiguous Q/E):", species.Z_count, file=obj)

    return # attributes added to species should be there after return


# print program name and version
print('===================================================================')
print(' program check_fasta_dir_walk.py, v1.0.0, Phil Wilmarth, OHSU 2020 ')
print('===================================================================')

# select a root folder
root_path = fasta_lib.get_folder(os.getcwd(), 'Select a Root folder')
if not root_path:
    sys.exit()     # cancel button repsonse

# create a log file to mirror screen output
_folder = root_path
log_obj = open(os.path.join(_folder, 'fasta_utilities.log'), 'a')
write = [None, log_obj]
fasta_lib.time_stamp_logfile('\n>>> starting: check_fasta.py', log_obj)

species_list = []
# process the FASTA files
for root, dirs, files in os.walk(root_path):
    for file in files:
        if file.endswith(".all.fa.gz"):
            species = Species()
            species.fasta_file = os.path.join(root,file)
            sub_folder = os.path.split(root)[-1]
            parts = sub_folder.split('__')
            species.common_name = parts[0]
            species.latin_name = parts[1]
            species.taxon = parts[2]
            species.release = file.split('__')[0]
            species.db_name = file.split('__')[2]
            try:
                fasta_checker(species, write)
                species_list.append(species)
            except FileNotFoundError:
                pass
            for obj in write:
                print(file=obj)

# finish up the log file
fasta_lib.time_stamp_logfile('>>> ending: check_fasta.py', log_obj)
log_obj.close()

# write the summary table
with open(os.path.join(root_path, 'database_analysis.txt'), 'wt') as fout:
    print(species_list[0].make_header(), file=fout)
    for species in species_list:
        print(species.make_row(), file=fout)

# end
