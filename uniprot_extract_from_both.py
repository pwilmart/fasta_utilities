"""'uniprot_extract_from_both.py' Written by Phil Wilmarth, OHSU.

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
# updtaed for Python 3 -PW 7/6/2017
import os
import sys
import fasta_lib_Py3 as fasta_lib

# set minimum sequence counts here
MIN_SEQUENCE_COUNT = 10      # minimum number of proteins per species
EXPAND_GROUPS = True        # controls expanding taxonomy groups (nodes)
MIN_GROUP_SEQ_COUNT = 50    # minimum if expanding a taxon group

# set this to True to simplify accession or False to keep unchanged
# NOTE: Cleaning accession/descriptions could cause some loss of information
CLEAN_ACCESSIONS = False
VERBOSE = True
MISMATCHES = False  # reports discrepancies between "spec_list.txt" and NCBI

# list species to extract by taxonomy number and name to use in filenames
taxon_dict = { 9606:'human',
               10090:'mouse',
               10116:'rat',
               559292:'yeast',
               419947:'M_tuberculosis_H37Ra',
               210007:'S_mutans_UA159',
               5811:'Toxoplasma_gondi',
               5141:'Neurospora_crassa',
               410289:'M_bovis_BCG_Pasteur_1173P2',
               246196:'M.smegmatis',
               9544:'Macaca.mulatta',
               9031:'chicken',
               5759:'E.histolytica',
               9615:'dog',
               3208:'Bryophyta',
               243243:'Mycobacterium_Avium',
               145481:'Physcomitrella_patens',
               5141:'Neurospora_Crassa',
               5693:'T.cruzi',
               5661:'L.donovani',
               7108:'Spodoptera.frugiperda',
               7227:'Drosophila.melanogaster',
               7091:'Bombyx.mori',
               5702:'T.brucei',
               9940:'sheep',
               8355:'Xenopus.laevis',
               161537: 'Bacillus.sp_pl-12',
               224308:'Bacillus.subtilis-168',
               10359:'Human_Cytomegalovirus',
               103930:'Rhesus_Cytomegalovirus',
               9823:'pig'}
##taxon_dict = { 145481:'Physcomitrella_patens'}


def main(taxon_dict):
    """Extracts entries by taxon ID from both Sprot and Trembl databases.
    """
    print('=============================================================================')
    print(' uniprot_extract_from_both.py, v.1.1.0, written by Phil Wilmarth, OHSU, 2017 ')
    print('=============================================================================')
    
    # get the UniProt folder and then get the sprot and trembl database names
    DB = []
    default = r'C:\Xcalibur\database'
    if not os.path.exists(default):
        default = os.getcwd()
    uniprot_folder = fasta_lib.get_folder(default, title_string='Select a UniProt download folder')
    if uniprot_folder == '': sys.exit() # cancel button response
    
    version = uniprot_folder.split('_')[-1]
    uniprot_db = 'uniprot'
    for files in os.listdir(uniprot_folder):
        if files.startswith('uniprot_') and files.endswith('.gz'):
            DB.append(os.path.join(uniprot_folder, files))
    if len(DB) != 2:
        print('WARNING: either sprot or trembl DB was missing')                      

    # create a log file to mirror screen output
    log_obj = open(os.path.join(uniprot_folder, 'fasta_utilities.log'), 'a')
    write = [None, log_obj]
    fasta_lib.time_stamp_logfile('\n>>> starting: uniprot_extract_from_both.py', log_obj)
    
    # make the smaller uniprot dictionaries
    (sci_to_taxon, id_to_taxon) = fasta_lib.make_uniprot_to_taxon(uniprot_folder)
    
    # make the more complete dictionary
    name_to_taxon = fasta_lib.make_all_names_to_taxon(uniprot_folder)
    
    # print the list of taxon numbers that will be extracted
    # NOTE: Any taxon numbers present in analysis text file will not be expanded.
    taxon_list = list(taxon_dict.items())
    taxon_list.sort()
    for obj in write:
        print('...extracting these taxon numbers:', file=obj)
        for i, t in enumerate(taxon_list):
            print('......(%s) taxon %s to file tagged with "%s"' % (i+1, t[0], t[1]), file=obj)
    
    # expand any group taxon numbers
    if EXPAND_GROUPS:
        fasta_lib.expand_species(uniprot_folder, 'uniprot', taxon_dict,
                                 MIN_SEQUENCE_COUNT, MIN_GROUP_SEQ_COUNT)
    
    # inititalize dictionaries and counters
    taxon_files, taxon_count, name_count = {}, {}, {}
    for taxon, name in taxon_dict.items():
        fname = uniprot_db+'_'+version+'_'+name+'.fasta'
        fname = os.path.join(uniprot_folder, fname)
        taxon_files[name] = fname
        taxon_count[taxon] = 0
        name_count[name] = 0
    
    # open the output filenames
    for name in taxon_files.keys():
        taxon_files[name] = open(taxon_files[name], 'w')
    
    # want to count extracted sequences from each database
    name_counter = {}
    number_counter = {}
    
    # loop over both databases and extract species
    duplicates = {}
    for i in range(len(DB)):
        prot_read = 0
        not_found = 0
        for value in taxon_dict.values():
            name_counter[value] = 0
        for key in taxon_dict.keys():
            number_counter[key] = 0
        
        # create a FastaReader object, initialize counters, and start reading
        uniprot_file = DB[i]
        x = fasta_lib.FastaReader(uniprot_file)
        prot = fasta_lib.Protein()
        for obj in write:
            print('...reading %s and extracting entries...' % (os.path.split(uniprot_file)[1],), file=obj)
        
        # NOTE: checking for errors will slow program execution, use if needed
        while x.readNextProtein(prot, check_for_errs=False):
            prot_read += 1
            if (prot_read % 500000) == 0:
                print('......(%s proteins read...)' % ("{0:,d}".format(prot_read),))
            (spec_id, spec_name) = fasta_lib.uniprot_parse_line(prot.accession + ' ' + prot.description)
            taxon = sci_to_taxon.get(spec_name, 0) # first choice mapping
            taxon2 = name_to_taxon.get(spec_name, 0) # alternative mapping
            if taxon == 0: # first choice not present
                if taxon2 == 0:
                    not_found += 1
                else:
                    taxon = taxon2 # use second choice
            else:
                if (taxon != taxon2) and (taxon2 > 0): # keep track of multiple taxon numbers
                    duplicates[spec_name] = (taxon, taxon2)
            if taxon_dict.get(taxon, False):
                if CLEAN_ACCESSIONS:
                    prot.parseUniProt()

                # taxon number matches, so write the protein to respective output file(s)
                name = taxon_dict[taxon]
                name_counter[name] += 1
                name_count[name] += 1
                taxon_count[taxon] += 1
                number_counter[taxon] += 1
                f = taxon_files[name]
                prot.printProtein(f)
        
        # print extraction stats for each database
        for obj in write:
            print('...%s protein entries in %s' %
                  ("{0:,d}".format(prot_read), os.path.split(DB[0])[1]), file=obj)
            print('...%s proteins had unknown taxon numbers' %
                  ("{0:,d}".format(not_found),), file=obj)
            numbers = list(number_counter.keys())
            numbers.sort()
            if VERBOSE:
                for j, number in enumerate(numbers):
                    if number_counter[number] > 0:
                        print('......(%s) taxon %s had %s proteins' %
                              (j+1, number, "{0:,d}".format(number_counter[number])), file=obj)
            names = list(name_counter.keys())
            names.sort()
            db_name = os.path.split(DB[i])[1]
            for j, name in enumerate(names):
                print('......(%s) %s %s proteins extracted' %
                      (j+1, "{0:,d}".format(name_counter[name]), name), file=obj)

    # close the extracted database files
    for f in taxon_files.values():
        f.close()
    
    # print list of mis-matched taxon number warnings
    if MISMATCHES:
        for i, (name, pair) in enumerate(duplicates.items()):
            for obj in write:
                print('......(%s) WARNING: %s and %s map to "%s"' % (i+1, pair[0], pair[1], name), file=obj)
    
    # print out the final summary stuff
    for obj in write:
        if VERBOSE:
            print('...combined taxon counts...', file=obj)
            numbers = list(taxon_count.keys())
            numbers.sort()
            for i, number in enumerate(numbers):
                if taxon_count[number] > 0:
                    print('......(%s) taxon %s had %s proteins' %
                          (i+1, number, "{0:,d}".format(taxon_count[number])), file=obj)
        print('...combined output file counts...', file=obj)
        for i, name in enumerate(names):
            print('......(%s) %s total proteins written to %s' %
                  (i+1, "{0:,d}".format(name_count[name]),
                   uniprot_db+'_'+version+'_'+name+'.fasta'), file=obj)
    
    fasta_lib.time_stamp_logfile('>>> ending: uniprot_extract_from_both.py', log_obj)
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


         
