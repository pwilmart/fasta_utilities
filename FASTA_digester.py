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
import fasta_lib_Py3 as fasta_lib


def fasta_digester(fasta_file, enzyme='trypsin', log=[None]):
    """Trypsin digests entries in a FASTA protein database.
        Call with FASTA filename, returns list of proteins with
        theoretical tryptic digest peptide lists
        Checks for duplicate accessions and (optional) valid characters.
    """
    print('==================================================================')
    print(' fasta_digester.py, v 1.1.3, written by Phil Wilmarth, OHSU, 2017 ')
    print('==================================================================')

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
    head = 0
    all_peptides = {}
    print('starting file reading:', time.ctime())
    
    # read proteins until EOF; NOTE: checking for errors slows program by factor of 3 or 4
    while f.readNextProtein(p, check_for_errs=False):

        # digest protein sequence (regex expression, low mass cutoff, high mas cutoff, 
        # minimum peptide length, maximum number of missed cleavages, type of masses)
        p.enzymaticDigest(regex, 500.0, 5000.0, 7, 2, 'mono')

        for pep in p.peptides:
            
            # mask I and L residues
            mass_spec_seq = re.sub(r'[IL]', 'j', pep.seq)
            
            # make dictionary of sequences and counts
            if all_peptides.get(mass_spec_seq):
                all_peptides[mass_spec_seq] += 1
            else:
                all_peptides[mass_spec_seq] = 1

        # count protein sequences
        prot += 1   
        if (prot % 500000) == 0:
            print('......(%s proteins read...)' % (prot,))
        
        # count number of header elements
        control_A = p.description.count(chr(1))
        head = head + control_A + 1
    
    # print number of proteins/headers and return peptide dictionary
    for obj in log:
        print('There are %s proteins in %s' %
              ("{0:,d}".format(prot), os.path.basename(fasta_file)), file=obj)
        if head > prot:
            print('There were %s header lines' % (head,), file=obj)

    return all_peptides
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
    log_obj = open(os.path.join(where, 'FASTA_digester_log.txt'), 'a')
    obj_list = [None, log_obj]
    
    # call digester function for each database
    for fasta_file in fasta_file_list:
        all_peptides = fasta_digester(fasta_file, log=obj_list)
        
        # compute unique and shared counts
        unique = [all_peptides[x] for x in all_peptides if all_peptides[x] == 1]
        shared = [all_peptides[x] for x in all_peptides if all_peptides[x] > 1]

        # find sequence with highest copy number
        sorted_peptides = sorted(all_peptides.items(), key=operator.itemgetter(1), reverse=True)
        
        # print some stats to the console
        ashared = numpy.array(shared)
        for obj in obj_list:
            print('\nDatabase:', os.path.basename(fasta_file), file=obj)
            print('Redundant counts:', file=obj)
            print('   total all peptides:', "{0:,d}".format(sum(all_peptides.values())), file=obj)
            print('   total shared peptides:', "{0:,d}".format(sum(shared)), file=obj)
            print('   total unique peptides:', "{0:,d}".format(sum(unique)), file=obj)
            print('   fraction unique:', round(100*float(sum(unique))/float(sum(all_peptides.values())),2), file=obj)
            print('   fraction shared:', round(100*float(sum(shared))/float(sum(all_peptides.values())),2), file=obj)
            print('Non-redundant counts:', file=obj)
            print('   total all non-redundant peptides:', "{0:,d}".format(len(all_peptides)), file=obj)
            print('   total shared peptides:', "{0:,d}".format(len(shared)), file=obj)
            print('   total unique peptides:', "{0:,d}".format(len(unique)), file=obj)
            print('   fraction unique:', round(100*float(len(unique))/float(len(all_peptides)),2), file=obj)
            print('   fraction shared:', round(100*float(len(shared))/float(len(all_peptides)),2), file=obj)
            print('Some shared peptide data:', file=obj)
            print('   max copy number:', "{0:,d}".format(max(shared)), sorted_peptides[0], file=obj)
            print('   min copy number:', min(shared), file=obj)
            print('   median copy number:', numpy.median(ashared), file=obj)
            print('   average copy number:', round(numpy.average(ashared), 2), file=obj)
            try:
                print('   upper quartile number:', round(numpy.percentile(ashared, 75, interpolation='higher'), 2), file=obj)
                print('   lower quartile number:', round(numpy.percentile(ashared, 25, interpolation='lower'), 2), file=obj)
            except TypeError:
                print('   upper quartile number:', round(numpy.percentile(ashared, 75), 2), file=obj)
                print('   lower quartile number:', round(numpy.percentile(ashared, 25), 2), file=obj)
            print(file=obj)

    log_obj.close()
    print('completed:', time.ctime())
                    
# end

