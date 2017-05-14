"""'uniprot_extract_from_one.py' Written by Phil Wilmarth, OHSU.
Copyright 2009, Oregon Health & Science University.
All Rights Reserved.

Permission to use, copy, modify, and distribute any part of this program
for non-profit scientific research or educational use, without fee, and
without a written agreement, is hereby granted, provided that the above
copyright notice, and this license agreement appear in all copies.
Inquiries regarding use of this software in commercial products or for
commercial purposes should be directed to:

Technology & Research Collaborations, Oregon Health & Science University,
2525 SW 1st Ave, Suite 120, Portland, OR 97210
Ph: 503-494-8200, FAX: 503-494-4729, Email: techmgmt@ohsu.edu.

IN NO EVENT SHALL OREGON HEALTH & SCIENCE UNIVERSITY BE LIABLE TO ANY
PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES,
INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE.  THE
SOFTWARE IS PROVIDED "AS IS", AND OREGON HEALTH &SCIENCE UNIVERSITY HAS
NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, OR ENHANCEMENTS.
OREGON HEALTH & SCIENCE UNIVERSITY MAKES NO REPRESENTATIONS NOR EXTENDS
WARRANTIES OF ANY KIND, EITHER IMPLIED OR EXPRESS, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A
PARTICULAR PURPOSE, OR THAT THE USE OF THE SOFTWARE WILL NOT INFRINGE
ANY PATENT, TRADEMARK OR OTHER RIGHTS.
"""
#
# set minimum sequence counts here
#
MIN_SEQUENCE_COUNT = 10      # minimum number of proteins per species
EXPAND_GROUPS = True        # controls expanding taxonomy groups (nodes)
MIN_GROUP_SEQ_COUNT = 50    # minimum if expanding a taxon group
#
# set this to True to simplify accession or False to keep unchanged
# NOTE: some protein description information may be lost during "cleaning"
#
CLEAN_ACCESSIONS = False
VERBOSE = True
MISMATCHES = False  # reports conflicts between "spec_list.txt" and NCBI taxons
#
# list species to extract by taxonomy number and name to use in filenames
#
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
               5693:'T.cruzi',
               5702:'T.brucei',
               5661:'L.donovani',
               7091:'Bombyz.mori',
               9031:'chicken',
               419947:'M.tuberculosis',
               246196:'M.smegmatis',
               9913:'Bovine',
               5141:'Neurospora_crassa',
               7227:'Drosophila.melanogaster'}
##taxon_dict = { 10360:'Human-cytomegalovirus_AD169'}
##taxon_dict = { 9986:'rabbit'}
##taxon_dict = { 244589:'SSV1'}
##taxon_dict = { 83333:'Ecoli'}
##taxon_dict = { 11686:'HIV_type1_groupM_subtypeB_isolateBRU-LAI'}
##taxon_dict = { 419947:'M.tuberculosis',
##               246196:'M.smegmatis'}
##taxon_dict = { 419947:'M_tuberculosis_H37Ra'}
##taxon_dict = { 7160:'Asian_Tiger_Mosquito' }
##taxon_dict = { 1309:'Streptococcus.muants' }
##taxon_dict = { 9940:'sheep' }
##taxon_dict = { 83333:'ecoli' }
##taxon_dict = { 5693:'T.cruzi' }
##taxon_dict = { 5702:'T.brucei' }
##taxon_dict = { 5661:'L.donovani' }
##taxon_dict = { 7091:'Bombyz.mori' }
##taxon_dict = { 9031:'chicken' }
##taxon_dict = { 9615:'dog' }
##taxon_dict = { 5759:'E.histolytica' }
##taxon_dict = { 9544:'Macaca.mulatta' }
##taxon_dict = { 419947:'M.tuberculosis',
##               246196:'M.smegmatis' }
##taxon_dict = { 243243:'M.avium',
##               1767:'M.Intracellulare'}
##taxon_dict = { 998088:'A.veronii' }
##taxon_dict = { 139:'Borrelia.burgdorferi' }
##taxon_dict = { 188937:'Methanosarcina.acetivorans' }
##taxon_dict = { 9913:'Bovine' }
##taxon_dict = { 7227:'Drosophila.melanogaster' }
##taxon_dict = { 5141:'Neurospora_crassa'}
#
#
#====================
def main(taxon_dict):
#====================
    """Main program to extract entries by taxon ID from uniprot databases.
    Extraction is from a single downloaded Sprot or Trembl database.
    """
    import fasta_lib
    import os
    import sys
    #
    print '============================================================================'
    print ' uniprot_extract_from_one.py, v.1.01, written by Phil Wilmarth, OHSU, 2009.'
    print '============================================================================'
    #
    # set some file paths and names
    #
    default = r'C:\Xcalibur\database'
    if not os.path.exists(default):
        default = os.getcwd()
    uniprot_file = fasta_lib.get_file(default, [('Zipped files', '*.gz'), \
                                                ('Fasta files', '*.fasta')], \
                                      title_string=\
                                      'Select an Sprot or Trembl database')
    if uniprot_file == '' : sys.exit() # cancel button repsonse
    uniprot_folder, uniprot_name = os.path.split(uniprot_file)
    version = uniprot_name.split('_')[-1]
    version = version.replace('.fasta.gz', '')
    uniprot_db = uniprot_name.split('_')[1]
    #
    # create a log file to mirror screen output
    #
    log_obj = open(os.path.join(uniprot_folder, 'fasta_utilities.log'), 'a')
    write = [None, log_obj]
    fasta_lib.time_stamp_logfile('\n>>> starting: uniprot_extract_from_one.py',\
                                 log_obj)
    #
    # make the smaller uniprot dictionaries
    #
    (sci_to_taxon, id_to_taxon) = fasta_lib.make_uniprot_to_taxon(uniprot_folder)
    #
    # make the more complete dictionary
    #
    name_to_taxon = fasta_lib.make_all_names_to_taxon(uniprot_folder)
    #
    # print the list of taxon numbers that will be extracted
    #
    taxon_list = taxon_dict.items()
    taxon_list.sort()
    for obj in write:
        print >>obj, '...extracting these taxon numbers:'
        for i, t in enumerate(taxon_list):
            print >>obj, '......(%s) taxon %s to file tagged with "%s"' % \
                  (i+1, t[0], t[1])
    #
    # expand any group taxon numbers
    # NOTE: Any taxon numbers present in "uniprot_fasta_analyze.txt" will not
    # be expanded.
    #
    if EXPAND_GROUPS:
        fasta_lib.expand_species(uniprot_folder, uniprot_db, taxon_dict, \
                                 MIN_SEQUENCE_COUNT, MIN_GROUP_SEQ_COUNT)
    #
    # inititalize dictionaries and counters
    #
    taxon_files = {}
    taxon_count = {}
    name_count = {}
    for taxon, name in taxon_dict.items():
        fname = uniprot_db+'_'+version+'_'+name+'.fasta'
        fname = os.path.join(uniprot_folder, fname)
        taxon_files[name] = fname
        taxon_count[taxon] = 0
        name_count[name] = 0
    #
    # open the output filenames
    #
    for name in taxon_files.keys():
        taxon_files[name] = open(taxon_files[name], 'w')
    #
    # create a FastaReader object, initialize counters, and start reading
    #
    x = fasta_lib.FastaReader(uniprot_file)
    prot = fasta_lib.Protein()
    prot_read = 0
    not_found = 0
    duplicates = {}
    for obj in write:
        print >>obj, '...reading %s and extracting entries...' % (uniprot_name,)
    #
    # checking for errors in sequences slows program execution, use as needed
    #
    while x.readNextProtein(prot, check_for_errs=False):
        prot_read += 1
        if (prot_read % 500000) == 0:
            print '......(%s proteins read...)' % (prot_read,)
        (spec_id, spec_name) = fasta_lib.uniprot_parse_line(prot.accession \
                                                            +' '+prot.description)    
        taxon = sci_to_taxon.get(spec_name, 0)
        taxon2 = name_to_taxon.get(spec_name, 0)
        if taxon == 0:
            if taxon2 == 0:
                not_found += 1
            else:
                taxon = taxon2
        else:
            if (taxon != taxon2) and (taxon2 > 0):
                duplicates[spec_name] = (taxon, taxon2)
        if taxon_dict.get(taxon, False):
            if CLEAN_ACCESSIONS:
                prot.parseUniProt()
            name = taxon_dict[taxon]
            name_count[name] += 1
            taxon_count[taxon] += 1
            f = taxon_files[name]
            prot.printProtein(f)
    for f in taxon_files.values():
        f.close()
    #
    # print list of mis-matching taxon number warnings
    #
    if MISMATCHES:
        items = duplicates.items()
        i = 1
        for name, pair in items:
            for obj in write:
                print >>obj, '......(%s) WARNING: %s and %s map to "%s"' % \
                      (i, pair[0], pair[1], name)
            i += 1
    #
    # print out the summary stuff
    #
    for obj in write:
        print >>obj, '...%s protein entries in %s' % (prot_read, uniprot_name)
        print >>obj, '...%s proteins had unknown taxon numbers' % (not_found,)
        if VERBOSE:
            numbers = taxon_count.keys()
            numbers.sort()
            for i, number in enumerate(numbers):
                if taxon_count[number] > 0:
                    print >>obj, '......(%s) taxon %s had %s proteins' % \
                          (i+1, number, taxon_count[number])
        print >>obj, '...output file summaries...'
        names = taxon_files.keys()
        names.sort()
        for i, name in enumerate(names):
            print >>obj, '......(%s) %s proteins extracted and written to %s' % \
                  (i+1, name_count[name], uniprot_db+'_'+version+'_'+name+'.fasta')
    #
    fasta_lib.time_stamp_logfile('>>> ending: uniprot_extract_from_one.py', log_obj)
    log_obj.close()
    #
    return
#
# check for command line launch and see if any arguments passed
#
import sys
import os
import fasta_lib
if __name__ == '__main__':
    if len(sys.argv) > 1:
        arg_dict = fasta_lib.taxon_cmd_line_checker(sys.argv)
        if arg_dict:
            main(arg_dict)
        else:
            sys.exit()
    else:
        main(taxon_dict)
    #
    try:    # wait for user to end program if not running from IDLE
        # what __file__ is under XP IDLE: 'C:\\Python26\\Lib\\idlelib\\idle.pyw'
        if not __file__.endswith('idle.pyw'):
            raw_input('\n...hit any key to end program...')
    except NameError:
        pass
#
# end
#


         
