"""'uniprot_extract_from_both.py' Written by Phil Wilmarth, OHSU.
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
# NOTE: Cleaning accession/descriptions could cause some loss of information
#
CLEAN_ACCESSIONS = False
VERBOSE = True
MISMATCHES = False  # reports discrepancies between "spec_list.txt" and NCBI
#
# list species to extract by taxonomy number and name to use in filenames
#
##taxon_dict = { 11676:'HIV-1',
##               11709:'HIV-2'}
##taxon_dict = { 388919:'Streptococcus.sanguinis_SK36',
##               1035195:'Corynebacterium durum'}
##taxon_dict = { 9606:'human',
##               10090:'mouse',
##               10116:'rat',
##               559292:'yeast',
##               419947:'M_tuberculosis_H37Ra',
##               210007:'S_mutans_UA159',
##               5811:'Toxoplasma_gondi',
##               5141:'Neurospora_crassa',
##               410289:'M_bovis_BCG_Pasteur_1173P2',
##               246196:'M.smegmatis',
##               9544:'Macaca.mulatta',
##               9031:'chicken',
##               5759:'E.histolytica',
##               9615:'dog',
##               3208:'Bryophyta',
##               243243:'Mycobacterium_Avium',
##               145481:'Physcomitrella_patens',
##               5141:'Neurospora_Crassa',
##               5693:'T.cruzi',
##               5661:'L.donovani',
##               7108:'Spodoptera.frugiperda',
##               7227:'Drosophila.melanogaster',
##               7091:'Bombyx.mori',
##               5702:'T.brucei',
##               9940:'sheep',
##               8355:'Xenopus.laevis',
##               161537: 'Bacillus.sp_pl-12',
##               224308:'Bacillus.subtilis-168',
##               10359:'Human_Cytomegalovirus',
##               103930:'Rhesus_Cytomegalovirus',
##               9823:'pig'}
taxon_dict = { 145481:'Physcomitrella_patens'}
##taxon_dict = { 680187:'Campylopterus',
##               1507371:'Campylopterus',
##               1315733:'Campylopterus',
##               689205:'Campylopterus',
##               472786:'Campylopterus',
##               304605:'Campylopterus',
##               926063:'Campylopterus',
##               304606:'Campylopterus'}
##taxon_dict = { 419947:'M_tuberculosis_H37Ra'}
##taxon_dict = { 1347420:'A_veronii_Hm21'}
##taxon_dict = { 243243:'M_avium_104'}
##taxon_dict = { 9823:'pig'}
##taxon_dict = { 9986:'rabbit',
##               9615:'dog'}
##taxon_dict = { 210007:'S_mutans_UA159'}
##taxon_dict = { 5811:'Toxoplasma_gondi'}
##taxon_dict = { 273057:'S_solfataricus_P2'}
##taxon_dict = { 5141:'Neurospora_crassa'}
##taxon_dict = { 410289:'M_bovis_BCG_Pasteur_1173P2'}
##taxon_dict = { 419947:'M.tuberculosis',
##               246196:'M.smegmatis'}
##taxon_dict = { 8355:'laevis',
##               8364:'tropicalis' }
##taxon_dict = { 37296:'Human_herpesvirus_8' }
##taxon_dict = { 224326:'Borrelia.burgdorferi.B31' }
##taxon_dict = { 11090:'Yellow_fever_virus_D17' }
##taxon_dict = { 9534:'Chlorocebus.aethiops' }
##taxon_dict = { 60711:'Chlorocebus.sabaeus' }
##taxon_dict = { 9544:'Macaca.mulatta' }
##taxon_dict = { 419947:'M.tuberculosis',
##               246196:'M.smegmatis' }
##taxon_dict = { 188937:'Methanosarcina.acetivorans' }
##taxon_dict = { 3208:'Bryophyta' }
##taxon_dict = { 9031:'chicken' }
##taxon_dict = { 5759:'E.histolytica' }
##taxon_dict = { 9615:'dog' }
##taxon_dict = { 243243:'Mycobacterium_Avium' }
##taxon_dict = { 9544:'Macaca mulatta' }
##taxon_dict = { 243243:'M.avium',
##               1767:'M.Intracellulare' }
##taxon_dict = { 145481:'Physcomitrella_patens' }
##taxon_dict = { 5141:'Neurospora_Crassa' }
##taxon_dict = { 1769:'M_leprae', 561304:'M_leprae_Br4923' }
##taxon_dict = { 5693:'T.cruzi' }
##taxon_dict = { 998088:'A.veronii' }
##taxon_dict = { 5661:'L.donovani' }
##taxon_dict = { 7108:'Spodoptera.frugiperda' }
##taxon_dict = { 7227:'Drosophila.melanogaster' }
##taxon_dict = { 7091:'Bombyx.mori' }
##taxon_dict = { 5702:'T.brucei' }
##taxon_dict = { 7160:'Asian_Tiger_Mosquito' }
##taxon_dict = { 1309:'Streptococcus.muants' }
##taxon_dict = { 9940:'sheep' }
##taxon_dict = { 8355:'Xenopus.laevis' }
##taxon_dict = { 936156:'Bacillus.subtilis-BSn5',
##taxon_dict = { 161537: 'Bacillus.sp_pl-12',
##               224308:'Bacillus.subtilis-168'}
##taxon_dict = { 1426:'Geobacillus.thermoglucosidas',
##               224308:'Bacillus.subtilis-168'}
##taxon_dict = { 184922:'Giardia.lamblia-50803',
##               598745:'Giardia.lamblia-50581',
##               658858:'Giardia.lamblia-p15',
##               419947:'Mycobacterium.tuberculosis-H37RV' }
##taxon_dict = { 10359:'Human_Cytomegalovirus' }
##taxon_dict = { 103930:'Rhesus_Cytomegalovirus' }
##taxon_dict = { 9823:'pig' }
##taxon_dict = { 116150: 'moths_buttterflies',
##               7091:   'moths_buttterflies',
##               13037:  'moths_buttterflies',
##               66420:  'moths_buttterflies',
##               29058:  'moths_buttterflies',
##               76194:  'moths_buttterflies',
##               7130:   'moths_buttterflies',
##               51655:  'moths_buttterflies' }
#
#
#====================
def main(taxon_dict):
#====================
    """Extracts entries by taxon ID from both Sprot and Trembl databases.
    """
    import fasta_lib
    import os
    import sys
    #
    print '============================================================================='
    print ' uniprot_extract_from_both.py, v.1.01, written by Phil Wilmarth, OHSU, 2009.'
    print '============================================================================='
    #
    # get the UniProt folder and then get the sprot and trembl database names
    #
    DB = []
    default = r'C:\Xcalibur\database'
    if not os.path.exists(default):
        default = os.getcwd()
    uniprot_folder = fasta_lib.get_folder(default, title_string=\
                                          'Select a UniProt download folder')
    if uniprot_folder == '': sys.exit() # cancel button response
    version = uniprot_folder.split('_')[-1]
    uniprot_db = 'uniprot'
    for files in os.listdir(uniprot_folder):
        if files.startswith('uniprot_') and files.endswith('.gz'):
            DB.append(os.path.join(uniprot_folder, files))
    if len(DB) != 2:
        print 'WARNING: either sprot or trembl DB was missing'                      
    #
    # create a log file to mirror screen output
    #
    log_obj = open(os.path.join(uniprot_folder, 'fasta_utilities.log'), 'a')
    write = [None, log_obj]
    fasta_lib.time_stamp_logfile('\n>>> starting: uniprot_extract_from_both.py',\
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
    #
    if EXPAND_GROUPS:
        fasta_lib.expand_species(uniprot_folder, 'uniprot', taxon_dict, \
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
    # want to count extracted sequences from each database
    #
    name_counter = {}
    number_counter = {}
    #
    # loop over both databases and extract species
    #
    duplicates = {}
    for i in range(len(DB)):
        prot_read = 0
        not_found = 0
        for value in taxon_dict.values():
            name_counter[value] = 0
        for key in taxon_dict.keys():
            number_counter[key] = 0
        #
        # create a FastaReader object, initialize counters, and start reading
        #
        uniprot_file = DB[i]
        x = fasta_lib.FastaReader(uniprot_file)
        prot = fasta_lib.Protein()
        for obj in write:
            print >>obj, '...reading %s and extracting entries...' %\
                  (os.path.split(uniprot_file)[1],)
        #
        # NOTE: checking for errors will slow program execution, use if needed
        #
        while x.readNextProtein(prot, check_for_errs=False):
            prot_read += 1
            if (prot_read % 500000) == 0:
                print '......(%s proteins read...)' % (prot_read,)
            (spec_id, spec_name) = fasta_lib.uniprot_parse_line(prot.accession \
                                                                + ' ' + \
                                                                prot.description)
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
                name_counter[name] += 1
                name_count[name] += 1
                taxon_count[taxon] += 1
                number_counter[taxon] += 1
                f = taxon_files[name]
                prot.printProtein(f)
        #
        # print extraction stats for each database
        #
        for obj in write:
            print >>obj, '...%s protein entries in %s' % \
                  (prot_read, os.path.split(DB[0])[1])
            print >>obj, '...%s proteins had unknown taxon numbers' % (not_found,)
            numbers = number_counter.keys()
            numbers.sort()
            if VERBOSE:
                for j, number in enumerate(numbers):
                    if number_counter[number] > 0:
                        print >>obj, '......(%s) taxon %s had %s proteins' % \
                              (j+1, number, number_counter[number])
            names = name_counter.keys()
            names.sort()
            db_name = os.path.split(DB[i])[1]
            for j, name in enumerate(names):
                print >>obj, '......(%s) %s %s proteins extracted' % \
                      (j+1, name_counter[name], name)        
    for f in taxon_files.values():
        f.close()
    #
    # print list of mis-matched taxon number warnings
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
    # print out the final summary stuff
    #
    for obj in write:
        if VERBOSE:
            print >>obj, '...combined taxon counts...'
            numbers = taxon_count.keys()
            numbers.sort()
            for i, number in enumerate(numbers):
                if taxon_count[number] > 0:
                    print >>obj, '......(%s) taxon %s had %s proteins' % \
                          (i+1, number, taxon_count[number])
        print >>obj, '...combined output file counts...'
        for i, name in enumerate(names):
            print >>obj, '......(%s) %s total proteins written to %s' % \
                  (i+1, name_count[name], uniprot_db+'_'+version+'_'+name+'.fasta')
    #
    fasta_lib.time_stamp_logfile('>>> ending: uniprot_extract_from_both.py', \
                                   log_obj)
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


         
