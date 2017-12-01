"""'Ensembl_fixer.py' Written by Phil Wilmarth, OHSU.

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

Simple Ensembl (post v83) FASTA protein parsing/repair program.
Truncates at stop codons, flags other odd things.
Parses and reformats description strings.
Written by Phil Wilmarth, OHSU, 2016
#
checks for duplicate accessions, too
added some additional testing to make header line parsing more robust -PW 12/1/2017
""" 
import os
import sys
import copy
import re
import argparse
import fasta_lib

                       
def parse_ensembl_header_line(line, all_tags):
    """Parses new format Ensembl FASTA header lines."""
    original = line
    parsed = {}
    tags = [x for x in all_tags if x in line]
    if not tags:    # might get empty tags if DB already fixed
        return line
    
    while line:
        line = line.strip() # do some trimming
        current_tag = tags.pop()    # get the next tag (from end to beginning)
        line, current_value = line.split(current_tag) # split out the next element
        parsed[current_tag] = current_value # save element in tag dictionary
        
    # build the desired description string
    if 'description:' not in parsed:
        parsed['description:'] = 'NO DESCRIPTION'
    extra = ' ('
    if 'gene:' in parsed and parsed['gene:']:
        if parsed['gene:'].startswith('ENS'):
            extra += 'g:' + str(float(re.sub("\D", "", parsed['gene:'])))
        else:
            extra += 'g:' + parsed['gene:']
    if 'transcript:' in parsed and parsed['transcript:']:
        if extra != ' (':
            extra += ', '
        if parsed['transcript:'].startswith('ENS'):
            extra += 't:' + str(float(re.sub("\D", "", parsed['transcript:'])))
        else:
            extra += 't:' + parsed['transcript:']
    if 'gene_symbol:' in parsed and parsed['gene_symbol:']:
        symbol = parsed['gene_symbol:']
        if extra != ' (':
            extra += ', '
        extra += 'gs:' + parsed['gene_symbol:']
    if extra == ' (':
        extra =  ''
    else:
        extra += ')'
    
    return parsed['description:'] + extra

def main(fasta_file, up_one=False):
    """Processes one Ensembl fasta file - reformats description lines, checks things.
    up_one determines where the new file is written.
    """
    # create the new database name
    original_fasta_file = os.path.basename(fasta_file)
    new_fasta_file = original_fasta_file.replace('.fasta', '_fixed.fasta')
    if new_fasta_file == original_fasta_file:
        new_fasta_file = original_fasta_file.replace('.fa', '_fixed.fasta')
    if new_fasta_file == original_fasta_file:
        print('WARNING! creating new file name failed')
        print('...make sure database is not compressed')
        return False
    if new_fasta_file.endswith('.gz'):
        new_fasta_file = new_fasta_file[:-3]
    if up_one:
        folder_name = os.path.dirname(os.path.dirname(fasta_file))
    else:
        folder_name = os.path.dirname(fasta_file)
    new_fasta_file = os.path.join(folder_name, new_fasta_file)

    # initializations
    proteins = []
    accessions = {}
    p = fasta_lib.Protein()
    pcount = 0      # sequence count
    dup_count = 0   # duplicate accession count
    stop_count = 0  # "*"
    gap_count = 0   # "-"
    no_met = 0      # does not start with M
    X_count = 0     # unknow AA
    B_count = 0     # N or D
    Z_count = 0     # Q or E
    J_count = 0     # I or L
    U_count = 0     # selenocysteine
    
    # set up the list of possible tags in header lines
    # this should probably be generalized somehow...
    all_tags = ['pep:', 'pep scaffold:', 'pep genescaffold:', 'pep chromosome:', 'pep contig:',
                'pep reftig:', 'pep supercontig:', 'pep ultracontig:', 'pep group:',
                'gene:', 'transcript:', 'gene_biotype:',
                'transcript_biotype:', 'gene_symbol:', 'description:']

    # read the sequences into a list
    f = fasta_lib.FastaReader(fasta_file)
    while f.readNextProtein(p, check_for_errs=True):
        pcount += 1
        
        # check if accession already seen
        if p.accession in accessions:
            dup_count += 1
            accessions[p.accession] += 1
            print('...WARNING: skipping duplicate accession:', p.accession)
            continue
        else:
            accessions[p.accession] = 1
        
        # clean up the description string
        p.new_desc = parse_ensembl_header_line(p.description, all_tags)
        
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
            p.new_desc = p.new_desc + ' (has gaps)'
        if 'B' in p.sequence:
            B_count += 1
            p.new_desc = p.new_desc + ' (has B)'
        if 'Z' in p.sequence:
            Z_count += 1
            p.new_desc = p.new_desc + ' (has Z)'
        if 'J' in p.sequence:
            J_count += 1
            p.new_desc = p.new_desc + ' (has J)'
        if 'U' in p.sequence:
            U_count += 1
            p.new_desc = p.new_desc + ' (has U)'
        if 'X' in p.sequence:
            X_count += 1
            p.new_desc = p.new_desc + ' (has unknown X)'
        
        # save the protein in list
        proteins.append(copy.deepcopy(p))

    # open the new protein fasta file and write out the proteins
    fixcount = 0
    file_obj = open(new_fasta_file, 'w')
    for p in proteins:
        if len(p.sequence) > 0:
            p.printProtein(file_obj)
        else:
            print('   empty sequence (stop codon at start):', p.accession)
        fixcount += 1
    file_obj.close()

    # print(out the report of oddball characters
    print("   Ensembl database:", os.path.basename(fasta_file))
    print("   translations that do not start with Met:", no_met)
    print("   translations that have premature stop codons:", stop_count)
    print("   translations that contain gaps:", gap_count)
    print("   translations that contain X (unknowns):", X_count)
    print("   translations that contain B:", B_count)
    print("   translations that contain Z:", Z_count)
    print("   translations that contain J:", J_count)
    print("   translations that contain U:", U_count)
    print("   total number of input sequences was:", pcount)
    print("   total number of sequences written was:", fixcount)
    print("   number of duplicate accessions was:", dup_count)

    return new_fasta_file
    # end

# setup stuff: check for command line args, etc.
if __name__ == '__main__':
    # set up command line arguments
    parser = argparse.ArgumentParser(description='Checks Ensembl databases and fixes descriptions.')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s version 1.1.1')
    parser.add_argument('file', help='list of FASTA files to process', nargs='*')

    args = parser.parse_args()
        
    if len(sys.argv) > 1:   # options set from command line
        fasta_files = args.files
    else:   # options set to hardcoded defaults if interactive mode or no passed commands
        fasta_files = []
    
    # if no FASTA files, browse to database file(s)
    if not fasta_files:
        database = r'C:\Xcalibur\database'
        if not os.path.exists(database):
            database = os.getcwd()
        fasta_files = fasta_lib.get_files(database,
                                          [('FASTA files', '*.fasta'), ('Zipped files', '*.gz'),
                                           ('fa files', '*.fa'), ('All files', '*.*')],
                                          'Select a FASTA database')
        if not fasta_files: sys.exit() # cancel button response

    # print version info, etc. (here because of the loop)    
    print('=================================================================')
    print(' Ensembl_fixer.py, v 1.1.1, written by Phil Wilmarth, OHSU, 2017 ')
    print('=================================================================')
    
    # call main function
    os.chdir('.')   # set location to where script lives - the contaminants should be there
    for fasta_file in fasta_files:
        try:
            main(fasta_file)
        except IOError:   # FastaReader class raises exception if file not found
            print('...WARNING: %s not found' % fasta_file)
            pass

# end
