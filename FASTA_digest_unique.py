"""'FASTA_digester.py' Written by Phil Wilmarth, OHSU.

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
import time
import re
import operator
import numpy
import copy
import fasta_lib


def fasta_digester(fasta_file, enzyme='trypsin', low_mass=500.0, high_mass=5000.0,
                   min_length=7, missed_cleavages=2, mass_type='mono', log=None):
    """Trypsin digests entries in a FASTA protein database.
        Call with FASTA filename, returns list of proteins with
        theoretical tryptic digest peptide lists
        Checks for duplicate accessions and (optional) valid characters.
    """
    print('=======================================================================')
    print(' fasta_digest_unique.py, v 1.0, written by Phil Wilmarth, OHSU, 2021 ')
    print('=======================================================================')

    # compile the regex for desired digestion
    if enzyme.upper() == 'No_enzyme'.upper():
        regex = re.compile(r".")
    elif enzyme.upper() == 'trypsin'.upper(): # checked
        regex = re.compile(r".(?:(?<![KR](?!P)).)*")
    elif enzyme.upper() == 'trypsin-P'.upper(): # checked
        regex = re.compile(r".(?:(?<![KR]).)*")
    elif enzyme.upper() == 'Lys-C'.upper(): # checked
        regex = re.compile(r".(?:(?<![K](?!P)).)*")
    elif enzyme.upper() == 'Lys-C-P'.upper(): # checked
        regex = re.compile(r".(?:(?<![K]).)*")
    elif enzyme.upper() == 'Lys-N'.upper(): # checked
        regex = re.compile(r".(?:(?![K]).)*")
    elif enzyme.upper() == 'Arg-C'.upper(): # checked
        regex = re.compile(r".(?:(?<![R](?!P)).)*")
    elif enzyme.upper() == 'Asp-N'.upper(): # checked
        regex = re.compile(r".(?:(?![D]).)*")
    elif enzyme.upper() == 'CNBr'.upper(): # checked
        regex = re.compile(r".(?:(?<![M]).)*")
    elif enzyme.upper() == 'Glu-C'.upper(): # checked
        regex = re.compile(r".(?:(?<![DE](?!P)).)*")
    elif enzyme.upper() == 'PepsinA'.upper(): # checked
        regex = re.compile(r".(?:(?<![FL](?!P)).)*")
    elif enzyme.upper() == 'chymotrypsin'.upper(): # checked
        regex = re.compile(r".(?:(?<![FWYL](?!P)).)*")
    else:
        print('...WARNING: Enzyme:', enzyme, 'not recognized')
        regex = None

    # create instances of reader object and protein object, initialize counters
    f = fasta_lib.FastaReader(fasta_file)
    p = fasta_lib.Protein()
    prot = 0
    proteins = []
    all_peptides = {}
    print('starting file reading:', time.ctime())

    # read proteins until EOF; NOTE: checking for errors slows program by factor of 3 or 4
    while f.readNextProtein(p, check_for_errs=False):

        # digest protein sequence (regex expression, low mass cutoff, high mas cutoff,
        # minimum peptide length, maximum number of missed cleavages, type of masses)
        p.enzymaticDigest(regex, low_mass, high_mass,
                          min_length, missed_cleavages, mass_type)

        # save all proteins that are read
        proteins.append(copy.copy(p))

        # count protein sequences
        prot += 1
        if (prot % 500000) == 0:
            print('......(%s proteins read...)' % (prot,))

    # print number of proteins/headers
    print('There are %s proteins in %s' %
          ("{0:,d}".format(prot), os.path.basename(fasta_file)), file=log)
        
    # make shared/unique status dictionary
    for p in proteins:
        for pep in p.peptides:
            # mask I and L residues
            mass_spec_seq = re.sub(r'[IL]', 'j', pep.seq)

            # make dictionary of sequences and counts
            if all_peptides.get(mass_spec_seq):
                all_peptides[mass_spec_seq].append(p.accession)
            else:
                all_peptides[mass_spec_seq] = [p.accession]

    keys = list(all_peptides.keys())
    print(keys[0], all_peptides[keys[0]])

    # print table (peptides from each protein, start, end, unique or not, protein list)
    print('\nAccession\tPeptide\tStart\tEnd\tMass\tMissed_Cleavages\tUnique\tOther_Proteins', file=log)

    for p in proteins:
        for pep in p.peptides:
            out_list = [p.accession]
            out_list += [pep.seq, str(pep.beg), str(pep.end), '%0.2f' % pep.mass, str(pep.missed)]
            mass_spec_seq = re.sub(r'[IL]', 'j', pep.seq)
            if len(all_peptides[mass_spec_seq]) == 1:
                out_list.append('TRUE')
            else:
                out_list.append('FALSE')

            acc_list = all_peptides[mass_spec_seq]
            if len(acc_list) == 1:
                acc_list = [' ']
            else:
                acc_list.remove(p.accession) 
                
            out_list += ['; '.join(acc_list)]

            # print table rows
            print('\t'.join(out_list), file=log)
            
    return
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
        fasta_file_list = fasta_lib.get_files(database,
                                              [('FASTA files', '*.fasta'),
                                               ('Zipped FASTA files', '*.gz'),
                                               ('All files', '*.*')],
                                              'Select a FASTA database')
        if not fasta_file_list: sys.exit()     # cancel button repsonse

    # set up file for results logging
    where = os.path.split(fasta_file_list[0])[0]
    log_obj = open(os.path.join(where, 'FASTA_digest_unique_log.txt'), 'wt')

    # call digester function for each database
    for fasta_file in fasta_file_list:
        # change digestion parameters here
        all_peptides = fasta_digester(fasta_file, enzyme='trypsin', low_mass=500.0, high_mass=5000.0,
                                      min_length=7, missed_cleavages=2, mass_type='mono', log=log_obj)

    log_obj.close()
    print('completed:', time.ctime())

# end
