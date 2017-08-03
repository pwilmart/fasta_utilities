"""'add_extras_and_reverse.py' Written by Phil Wilmarth, OHSU.

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
# fixed bug in the sequence reversing - May 2012 PW
# now compatible with python 3.6 - July 2017 PW

import os
import sys
import fasta_lib_Py3 as fasta_lib

# flag to parse (clean) accessions and descriptions
# NOTE: cleaning accessions/descriptions may cause some loss of information.
CLEAN_ACCESSIONS = False

# if RefSeq accesions are preferred, set to True (False returns gi numbers) 
REF_SEQ_ONLY = True

# flag to keep IPI gene identifiers (True) or not (False)
KEEP_IPI_GENE_ID = False

# flag to keep UniProt identifier (True) or more-stable accession (False)
KEEP_UNIPROT_ID = True

# flags to make different output databases
MAKE_SEPARATE_FORWARD = False
MAKE_SEPARATE_REVERSED = False
MAKE_SEPARATE_BOTH = True


def fasta_add_extras(extra_file, fasta_file, output_file):
    """Adds contaminants and reverses entries in a FASTA protein database.
        Called with FASTA filename.  Reversed DB written to same location.
        Options for separate or concatenated output files.
    """
    decoy_string = 'REV_'   # the string to denote decoy sequences
    
    print('=========================================================================')
    print(' add_extras_and_reverse.py, v1.1.0, written by Phil Wilmarth, OHSU, 2017 ')
    print('=========================================================================')
    
    # open the "forward" and "reversed" output files
    _file = os.path.splitext(output_file)[0] + '.fasta'
    for_name = _file.replace('.fasta', '_for.fasta')
    for_file_obj = open(for_name, 'w')
    rev_name = _file.replace('.fasta', '_rev.fasta')
    rev_file_obj = open(rev_name, 'w')
    
    # create a log file to mirror screen output
    _folder = os.path.split(fasta_file)[0]
    log_obj = open(os.path.join(_folder, 'fasta_utilities.log'), 'a')
    write = [None, log_obj]
    fasta_lib.time_stamp_logfile('\n>>> starting: add_extras_and_reverse.py', log_obj)
    
    # create instances of reader object and protein object
    # add the extra sequences, accessions changed to "EXTRA_dddd"
    # NOTE: can only add up to 9999 sequences without code modification...
    prot = fasta_lib.Protein()
    pcount = 0
    f = fasta_lib.FastaReader(extra_file)
    
    # turn on error checking for extra sequences
    while f.readNextProtein(prot, check_for_errs=True):
        pcount += 1
        
        # try to clean up original accessions
        if CLEAN_ACCESSIONS:
            if prot.accession.startswith('gi|'):
                prot.parseNCBI(REF_SEQ_ONLY)
            elif prot.accession.startswith('sp|') or prot.accession.startswith('tr|'):
                prot.parseUniProt(KEEP_UNIPROT_ID)
            elif prot.accession.startswith('IPI:'):
                prot.parseIPI(KEEP_IPI_GENE_ID)
            else:
                pass
        
        # add old accession to description and make new accession
        # write original sequence to "forward" file and reversed to "reverse"
        prot.new_desc = '[%s] %s' % (prot.new_acc, prot.new_desc)
        prot.new_acc = 'EXTRA_%04d' % (pcount,)
        prot.accession = 'EXTRA_%04d' % (pcount,)
        prot.printProtein(for_file_obj)
        rev = prot.reverseProtein(decoy_string)
        rev.printProtein(rev_file_obj)
    for obj in write:
        print('...there were %s extra sequences in %s' % (pcount, os.path.split(extra_file)[1]), file=obj)            
    
    # now add the contaminants
    try:
        if os.path.exists('all_contams_fixed.fasta'):
            contams_file = 'all_contams_fixed.fasta'
        else:
            path = os.path.split(fasta_file)[0]
            contams_file = os.path.join(path, 'all_contams_fixed.fasta')
        f = fasta_lib.FastaReader(contams_file)
        contams = 0
        while f.readNextProtein(prot, check_for_errs=True):
            pcount += 1
            contams += 1
            if CLEAN_ACCESSIONS:
                prot.parseCONT()
            
            # write sequences to respective files
            prot.printProtein(for_file_obj)
            rev = prot.reverseProtein(decoy_string)
            rev.printProtein(rev_file_obj)
        for obj in write:
            print('...there were %s contaminant entries in %s' % (contams, contams_file), file=obj)            
    except:
        for obj in write:
            print('...WARNING: "all_contams_fixed.fasta" not found!', file=obj)
    
    # read proteins, clean up accessions, decriptions until EOF
    # write proteins to "forward" and "reversed" files
    f = fasta_lib.FastaReader(fasta_file)
    
    # checking for errors can slow program execution by factor of 3-4
    # Reading and writing sequences will always remove spaces and blank lines
    while f.readNextProtein(prot, check_for_errs=False):
        pcount += 1
        if CLEAN_ACCESSIONS:
            if prot.accession.startswith('gi|'):
                prot.parseNCBI(REF_SEQ_ONLY)
            elif prot.accession.startswith('sp|') or prot.accession.startswith('tr|'):
                prot.parseUniProt(KEEP_UNIPROT_ID)
            elif prot.accession.startswith('IPI:'):
                prot.parseIPI(KEEP_IPI_GENE_ID)
            else:
                pass
        
        prot.printProtein(for_file_obj)    # write to "forward" file
        rev = prot.reverseProtein(decoy_string)
        rev.printProtein(rev_file_obj)   # write to "reversed" file
    
    # make concatenated output file if desired and print summary stats
    if MAKE_SEPARATE_BOTH:
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
            print('...%s total proteins written to %s' % (2*pcount, os.path.split(both_name)[1]), file=obj)
    
    if MAKE_SEPARATE_FORWARD:
        for obj in write:
            print('...%s proteins written to %s' % (pcount, os.path.split(for_name)[1]), file=obj)
    if MAKE_SEPARATE_REVERSED:
        for obj in write:
            print('...%s proteins reversed and written to %s' % (pcount, os.path.split(rev_name)[1]), file=obj)
    
    # close files and delete unwanted files
    for_file_obj.close()
    rev_file_obj.close()
    fasta_lib.time_stamp_logfile('>>> ending: add_extras_and_reverse.py', log_obj)
    log_obj.close()
    if not MAKE_SEPARATE_FORWARD:
        os.remove(for_name)
    if not MAKE_SEPARATE_REVERSED:
        os.remove(rev_name)
    return


# setup stuff: check for command line args, etc.
if __name__ == '__main__':    
    extra_file = ''
    fasta_file = ''
    output_file = ''
    
    # check if database names were passed on command line
    if len(sys.argv) == 4:
        if os.path.exists(sys.argv[1]):
            extra_file = sys.argv[1]
        if os.path.exists(sys.argv[2]):
            fasta_file = sys.argv[2]
        if os.path.exists(sys.argv[3]):
            output_file = sys.argv[3]
    
    # if not, get database files with dialog boxes
    else:
        if len(sys.argv) > 1:
            for i, db in enumerate([extra_file, fasta_file, output_file]):
                if not db:                    
                    print('...WARNING: %s not found...' % (sys.argv[i+1],))
        database = r'C:\Xcalibur\database'
        if not os.path.exists(database):
            database = os.getcwd()
            
        print('Select the FASTA file with extra sequences')
        extra_file = fasta_lib.get_file(database,
                                        [('FASTA files', '*.fasta'), ('All files', '*.*')],
                                        'Select Extra Sequences (FASTA format)')
        if extra_file == '': sys.exit() # cancel button response
        
        extra_name = os.path.split(extra_file)[1]
        extra_name = extra_name.split('.fasta')[0]        
        print('Select the main FASTA file')
        fasta_file = fasta_lib.get_file(database, [('FASTA files', '*.fasta'),
                                                   ('GZipped files', '*.gz'),
                                                   ('All files', '*.*')],
                                        'Select FASTA database file')
        if fasta_file == '': sys.exit() # cancel button response
        
        default = os.path.split(fasta_file)[0]
        fasta_name = os.path.split(fasta_file)[1]
        default_file = extra_name + '_' + fasta_name
        print('Select location and name for new FASTA database')
        output_file = fasta_lib.save_file(default, [('FASTA files', '*.fasta'),
                                                    ('All files', '*.*')],
                                          default_file, 'Output filename and location')
        if output_file == '': sys.exit() # cancel button response        
    
    # call add extras and reverse function
    fasta_add_extras(extra_file, fasta_file, output_file)

# end
