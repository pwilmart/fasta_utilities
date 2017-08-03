"""'uniprot_extract_from_one.py' Written by Phil Wilmarth, OHSU.

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
import fasta_lib_Py3 as fasta_lib

# set minimum sequence counts here
MIN_SEQUENCE_COUNT = 10      # minimum number of proteins per species
EXPAND_GROUPS = True        # controls expanding taxonomy groups (nodes)
MIN_GROUP_SEQ_COUNT = 50    # minimum if expanding a taxon group

# set this to True to simplify accession or False to keep unchanged
# NOTE: some protein description information may be lost during "cleaning"
CLEAN_ACCESSIONS = False
VERBOSE = True
MISMATCHES = False  # reports conflicts between "spec_list.txt" and NCBI taxons

# list species to extract by taxonomy number and name to use in filenames
##taxon_dict = { 9606: 'human',
##               10090: 'mouse',
##               10116: 'rat',
##               559292: 'yeast',
##               83333: 'Ecoli',
##               9913: 'Bovine'}
taxon_dict = { 9606: 'human',
               10090: 'mouse',
               10116: 'rat',
               559292: 'yeast',
               83333: 'Ecoli',
               419947:'M.tuberculosis',
               246196:'M.smegmatis',
               1309:'Streptococcus.muants',
               9940:'sheep',
               7091:'Bombyz.mori',
               9031:'chicken',
               419947:'M.tuberculosis',
               246196:'M.smegmatis',
               9913:'Bovine',
               5141:'Neurospora_crassa',
               7227:'Drosophila.melanogaster'}


def main(taxon_dict):
    """Main program to extract entries by taxon ID from uniprot databases.
    Extraction is from a single downloaded Sprot or Trembl database.
    """
    print('============================================================================')
    print(' uniprot_extract_from_one.py, v.1.1.0, written by Phil Wilmarth, OHSU, 2017 ')
    print('============================================================================')
    
    # set some file paths and names
    default = r'C:\Xcalibur\database'
    if not os.path.exists(default):
        default = os.getcwd()
    uniprot_file = fasta_lib.get_file(default,
                                      [('Zipped files', '*.gz'), ('Fasta files', '*.fasta')],
                                      title_string = 'Select an Sprot or Trembl database')
    if uniprot_file == '' : sys.exit() # cancel button repsonse
    
    uniprot_folder, uniprot_name = os.path.split(uniprot_file)
    version = uniprot_name.split('_')[-1]
    version = version.replace('.fasta.gz', '')
    uniprot_db = uniprot_name.split('_')[1]
    
    # create a log file to mirror screen output
    log_obj = open(os.path.join(uniprot_folder, 'fasta_utilities.log'), 'a')
    write = [None, log_obj]
    fasta_lib.time_stamp_logfile('\n>>> starting: uniprot_extract_from_one.py', log_obj)
    
    # make the smaller uniprot dictionaries
    (sci_to_taxon, id_to_taxon) = fasta_lib.make_uniprot_to_taxon(uniprot_folder)
    
    # make the more complete dictionary
    name_to_taxon = fasta_lib.make_all_names_to_taxon(uniprot_folder)
    
    # print the list of taxon numbers that will be extracted
    taxon_list = list(taxon_dict.items())
    taxon_list.sort()
    for obj in write:
        print('...extracting these taxon numbers:', file=obj)
        for i, t in enumerate(taxon_list):
            print('......(%s) taxon %s to file tagged with "%s"' % (i+1, t[0], t[1]), file=obj)

    # expand any group taxon numbers
    # NOTE: Any taxon numbers present in analysis text file will not be expanded.
    if EXPAND_GROUPS:
        fasta_lib.expand_species(uniprot_folder, uniprot_db, taxon_dict, MIN_SEQUENCE_COUNT, MIN_GROUP_SEQ_COUNT)
    
    # inititalize dictionaries and counters
    taxon_files, taxon_count, name_count = {}, {}, {}
    for taxon, name in taxon_dict.items():
        fname = uniprot_db + '_' + version + '_' + name + '.fasta'
        fname = os.path.join(uniprot_folder, fname)
        taxon_files[name] = fname
        taxon_count[taxon] = 0
        name_count[name] = 0
    
    # open the output filenames
    for name in taxon_files.keys():
        taxon_files[name] = open(taxon_files[name], 'w')
    
    # create a FastaReader object, initialize counters, and start reading
    x = fasta_lib.FastaReader(uniprot_file)
    prot = fasta_lib.Protein()
    prot_read = 0
    not_found = 0
    duplicates = {}
    for obj in write:
        print('...reading %s and extracting entries...' % (uniprot_name,), file=obj)
    
    # checking for errors in sequences slows program execution, use as needed
    while x.readNextProtein(prot, check_for_errs=False):
        prot_read += 1
        if (prot_read % 500000) == 0:
            print('......(%s proteins read...)' % ("{0:,d}".format(prot_read),))
        (spec_id, spec_name) = fasta_lib.uniprot_parse_line(prot.accession + ' ' + prot.description)    
        taxon = sci_to_taxon.get(spec_name, 0) # first choice mapping
        taxon2 = name_to_taxon.get(spec_name, 0) # alternative mapping
        if taxon == 0:  # first choice not present
            if taxon2 == 0:
                not_found += 1
            else:
                taxon = taxon2 # use second choice
        else:
            if (taxon != taxon2) and (taxon2 > 0): #keep track of multiple taxon numbers
                duplicates[spec_name] = (taxon, taxon2)
        if taxon_dict.get(taxon, False):
            if CLEAN_ACCESSIONS:
                prot.parseUniProt()

            # taxon number matches, so write the protein to the respective file
            name = taxon_dict[taxon]
            name_count[name] += 1
            taxon_count[taxon] += 1
            f = taxon_files[name]
            prot.printProtein(f)

    # close the extracted database files
    for f in taxon_files.values():
        f.close()
    
    # print list of mis-matching taxon number warnings
    if MISMATCHES:
        for i, (name, pair) in enumerate(duplicates.items()):
            for obj in write:
                print('......(%s) WARNING: %s and %s map to "%s"' % (i+1, pair[0], pair[1], name), file=obj)
    
    # print out the summary stuff
    for obj in write:
        print('...%s protein entries in %s' % ("{0:,d}".format(prot_read), uniprot_name), file=obj)
        print('...%s proteins had unknown taxon numbers' % (not_found,), file=obj)
        if VERBOSE:
            numbers = list(taxon_count.keys())
            numbers.sort()
            for i, number in enumerate(numbers):
                if taxon_count[number] > 0:
                    print('......(%s) taxon %s had %s proteins' %
                          (i+1, number, "{0:,d}".format(taxon_count[number])), file=obj)
        print('...output file summaries...', file=obj)
        names = list(taxon_files.keys())
        names.sort()
        for i, name in enumerate(names):
            print('......(%s) %s proteins extracted and written to %s' %
                  (i+1, "{0:,d}".format(name_count[name]),
                   uniprot_db + '_' + version + '_' + name + '.fasta'), file=obj)
    
    fasta_lib.time_stamp_logfile('>>> ending: uniprot_extract_from_one.py', log_obj)
    log_obj.close()
    return

# check for command line launch and see if any arguments passed
if __name__ == '__main__':
    if len(sys.argv) > 1:
        arg_dict = fasta_lib.taxon_cmd_line_checker(sys.argv)
        if arg_dict:
            main(arg_dict)
        else:
            sys.exit()
    else:
        main(taxon_dict)

# end


         
