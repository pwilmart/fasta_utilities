"""'remove_duplicates.py' Written by Phil Wilmarth, OHSU.

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
import copy
import fasta_lib_Py3 as fasta_lib


def find_identities(prot, candidates, skip, fasta_file, out_obj):
    """Test sequences in candiates list for identity to prot.
    """
    ident = 0
    acc = prot.accession
    seq = prot.sequence
    for p in candidates[acc]:
        if seq == p.sequence:
            ident += 1
            skip[p.accession] = True
            prot.new_desc = prot.new_desc + chr(1) + p.accession + ' ' + p.description
            if ident == 1:
                print('\n...%s "%s" same as:' % (acc, prot.description[:60]), file=out_obj)
                print('......%s "%s"' % (p.accession, p.description[:60]), file=out_obj)
            else:
                print('......%s "%s"' % (p.accession, p.description[:60]), file=out_obj)
    return ident

def main(fasta_file):
    """Checks entries in a FASTA protein database for identical duplicates.
        Call with FASTA filename, returns a couple of dictionaries
    """
    print('====================================================================')
    print(' remove_duplicates.py, v1.1.0, written by Phil Wilmarth, OHSU, 2017 ')
    print('====================================================================')
    print('\n   Analysis results will be written to "duplicates.txt"')
    
    # set up the output file names, etc.
    folder = os.path.split(fasta_file)[0]
    out_file = os.path.join(folder, 'duplicates.txt')
    out_obj = open(out_file, 'a')
    for end in ['.fasta', '.fasta.gz', '.fa.gz']:
        if fasta_file.endswith(end):
            nr_database = fasta_file.replace(end, '_nonredun.fasta')
    if (not nr_database) or (nr_database == fasta_file):
        nr_database = fasta_file + '_nonredun.fasta'
    nr_obj = open(nr_database, 'w')
    write = [None, out_obj]
    fasta_lib.time_stamp_logfile('\n>>> starting: check_for_duplicates.py', out_obj)
    #
    
    # create instances of reader object and protein object, initialize counters
    f = fasta_lib.FastaReader(fasta_file)
    p = fasta_lib.Protein()
    prot, head, dup = 0, 0, 0
    candidates = {}     # dictionary of acc:protein lists
    conflicts = {}      # keeps track of seq len and MW
    
    # read proteins until EOF
    while f.readNextProtein(p, check_for_errs=False):
        prot += 1
        control_A = p.description.count(chr(1))
        head = head + control_A + 1
        dup_data = (p.seqlenProtein(), p.molwtProtein())
        duplicate = conflicts.get(dup_data, False)
        if duplicate:
            dup += 1
            if candidates.get(duplicate, False):
                candidates[duplicate].append(copy.deepcopy(p))
            else:
                candidates[duplicate] = [copy.deepcopy(p)]
        else:
            conflicts[dup_data] = p.accession
    
    # get list of proteins to test for identity
    to_test = {}
    for key in candidates.keys():
        to_test[key] = True    
    
    # copy proteins to "nr" file, checking for duplicates
    dup = 0
    skip = {}
    f = fasta_lib.FastaReader(fasta_file)
    p = fasta_lib.Protein()
    print('Processing:', fasta_file, file=out_obj)    # header line to log file
    while f.readNextProtein(p, check_for_errs=False):
        if to_test.get(p.accession, False):
            dup += find_identities(p, candidates, skip, fasta_file, out_obj)
        if skip.get(p.accession, False):
            continue
        p.printProtein(nr_obj)
    
    for obj in [None, out_obj]:
        print('\nThere were', prot, 'total sequences in:', os.path.basename(fasta_file), file=obj)
        print('There were', dup, 'identical sequences removed\n\n', file=obj)
        try:
            obj.close()
        except AttributeError:
            pass
    return


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
    main(fasta_file)

# end
