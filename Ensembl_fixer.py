
"""Simple Ensembl (post v83) FASTA protein parsing/repair program.
Truncates at stop codons, flags other odd things.
Parses and reformats description strings.
Written by Phil Wilmarth, OHSU, 2016
#
checks for duplicate accessions, too
""" 

import os
import copy
import fasta_lib

# set up the list of possible tags in header lines
all_tags = all_tags = ['pep:', 'pep scaffold:', 'pep chromosome:', 'gene:', 'transcript:', 
                       'gene_biotype:', 'transcript_biotype:', 'gene_symbol:', 'description:']
                       
def parse_ensembl_header_line(line):
    """Parses new format Ensembl FASTA header lines."""
    parsed = {}
    tags = [x for x in all_tags if x in line]
    while line:
        line = line.strip() # do some trimming
        current_tag = tags.pop()    # get the next tag (from end to beginning)
        line, current_value = line.split(current_tag) # split out the next element
        parsed[current_tag] = current_value # save element in tag dictionary
        
    # build the desired description string
    if 'description:' not in parsed:
        parsed['description:'] = 'NO DESCRIPTION'
    extra = ' ('
    if 'gene:' in parsed:
        extra += 'g:' + str(float(parsed['gene:'][7:]))
    if 'transcript:' in parsed:
        if extra != ' (':
            extra += ', '
        extra += 't:' + str(float(parsed['transcript:'][7:]))
    if 'gene_symbol:' in parsed:
        symbol = parsed['gene_symbol:']
        if extra != ' (':
            extra += ', '
        extra += 'gene:' + parsed['gene_symbol:']
    if extra == ' (':
        extra =  ''
    else:
        extra += ')'
    
    return parsed['description:'] + extra
    
# print program name and version
print('==========================================================')
print(' program Ensembl_fixer.py, v1.0, Phil Wilmarth, OHSU 2016 ')
print('==========================================================')

# browse to the database
database_path = r"C:\Google_Drive\Databases"
if not os.path.exists(database_path):
    database_path = os.getcwd()
fasta_file = fasta_lib.get_file(database_path, [('FASTA files', '*.fasta')],
                                'Select an Ensembl FASTA database')
if fasta_file == '': sys.exit()     # cancel button repsonse

# create the new database name
new_fasta_file = os.path.basename(fasta_file)
new_fasta_file = new_fasta_file.replace('.fasta', '_fixed.fasta')
new_fasta_file = os.path.join(os.path.dirname(fasta_file), new_fasta_file)

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
    p.new_desc = parse_ensembl_header_line(p.description)
    
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
        print('zero sequence:', p.accession)
    fixcount += 1
file_obj.close()

# print out the report of oddball characters
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
#
# end
#
