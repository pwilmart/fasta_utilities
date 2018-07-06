"""'check_for_duplicates.py' Written by Phil Wilmarth, OHSU.

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
# updated for Python 3 -PW 7/6/2017

import os
import sys
import copy
import fasta_lib

def main(fasta_file):
    """Checks entries in a FASTA protein database for identical duplicates.
        Call with FASTA filename, returns a couple of dictionaries
    """
    print('=======================================================================')
    print(' check_for_duplicates.py, v1.1.0, written by Phil Wilmarth, OHSU, 2017 ')
    print('=======================================================================')

    # set up the output file names, etc.
    folder = os.path.split(fasta_file)[0]
    out_file = os.path.join(folder, 'duplicates.txt')
    out_obj = open(out_file, 'w')

    # create instances of reader object and protein object, initialize counter
    f = fasta_lib.FastaReader(fasta_file)
    p = fasta_lib.Protein()
    prot = 0
    head = 0
    dup = 0
    conflict = {}

    # read proteins until EOF
    while f.readNextProtein(p, check_for_errs=False):
        prot += 1
        control_A = p.description.count(chr(1))
        head = head + control_A + 1
        dup_data = (p.seqlenProtein(), p.molwtProtein())
        value = p.accession
        duplicate = conflict.get(dup_data, False)
        if duplicate:
            dup += 1
            print('...WARNING: (protein no. %s) %s may be same as %s' % (prot, p.accession, duplicate), file=out_obj)
        else:
            conflict[dup_data] = value

    # print result and return
    print('...there are %s proteins in %s' % (prot, os.path.split(fasta_file)[1]))
    if head > prot:
        print('...there were %s header lines...' % (head,))
    print('...there were %s possible duplicates...' % (dup,))
    out_obj.close()

    # rewind out file and build new dictionaries
    out_obj = open(out_file, 'r')
    conflict = {}
    to_save = {}
    while True:
        line = out_obj.readline()
        if not line:
            break
        else:
            line = line.split()
        conflict[line[4]] = line[9]
        to_save[line[4]] = True
        to_save[line[9]] = True
    out_obj.close()

    # read in and save the proteins that might be duplicates of each other
    candidates = []
    i = 0
    dup2 = 0
    index = {}
    f = fasta_lib.FastaReader(fasta_file)
    p = fasta_lib.Protein()
    while f.readNextProtein(p, check_for_errs=False):
        if to_save.get(p.accession, False):
            candidates.append(copy.deepcopy(p))
            index[p.accession] = i
            i += 1
    if len(candidates) == 0:
        print('to_save_dictionary:', to_save)
        print('bailing out in middle')
        return(candidates, index)

    # look deeper to see if candidates are actually duplicates
    out_obj = open(out_file, 'a')
    print('\n========================================\n', file=out_obj)
    exact_dup = 0
    accessions = list(conflict.keys())
    accessions.sort()
    for acc in accessions:
        p_dup = candidates[index[acc]]
        dup_acc = conflict[acc]
        p_ref = candidates[index[dup_acc]]
        if p_ref.sequence == p_dup.sequence:
            exact_dup += 1
            print('...(%s) WARNING: %s exact match to %s' % (exact_dup, p_ref.accession, p_dup.accession), file=out_obj)
    print('...number of exact matches was', exact_dup)

    return(candidates, index)
    # end


# setup stuff: check for command line args, etc.
if __name__ == '__main__':

    # check if database name passed on command line
    if len(sys.argv) > 1 and os.path.exists(sys.argv[1]):
        fasta_file = sys.argv[1]

    # if not, browse to database file
    else:
        if len(sys.argv) > 1:
            print('...WARNING: %s not found...' % (sys.argv[1],))
        database = r'C:\Xcalibur\database'
        if not os.path.exists(database):
            database = os.getcwd()
        fasta_file = fasta_lib.get_file(database,
                                        [('FASTA files', '*.fasta'), ('Zipped FASTA files', '*.gz'), ('All files', '*.*')],
                                        'Select a FASTA database')
        if fasta_file == '': sys.exit()     # cancel button repsonse

    # call main function
    candidates, index = main(fasta_file)

# end
