"""'uniprot_get_analyze.py' Written by Phil Wilmarth, OHSU.
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
# minimum sequence count cutoff for output table -
# there are MANY species in Trembl with 5 or fewer protein sequences
#
min_sequence_count = [0, 10]     # [sprot, trembl]
#
#
#==============================
def main(DB, folder, versions):
#==============================
    """Analyzes the species names in both UniProt databases.

    Arguments:
    "folder" is full path name to UniProt DBs.  "versions" is
    a dictionary of version numbers for file naming.  No return values.

    Saves a summary text file that can be loaded into EXCEL.
    """
    import os
    import fasta_lib
    #
    print '======================================================================'
    print ' uniprot_get_analyze.py, v1.01, written by Phil Wilmarth, OHSU, 2009.'
    print '======================================================================'
    #
    # create a log file to mirror screen output
    #
    log_obj = open(os.path.join(folder, 'fasta_utilities.log'), 'a')
    write = [None, log_obj]
    fasta_lib.time_stamp_logfile('\n>>> starting: uniprot_get_analyze.py', \
                                 log_obj)
    #
    # make sure the files are present or download if not
    #
    for db in DB:
        fasta_lib.download_uniprot(db, folder, versions)
    #
    # make dictionaries of species names (or accession IDs) to taxon IDs
    #
    (sci_to_taxon, id_to_taxon) = fasta_lib.make_uniprot_to_taxon(folder)
    #
    # get more complete list of names to taxon number from ncbi data
    #
    name_to_taxon = fasta_lib.make_all_names_to_taxon(folder)
    #
    # make species frequency dictionary
    #
    for i in range(len(DB)):
        fname = 'uniprot_%s_%s.fasta.gz' % (DB[i], versions[DB[i]],)
        db_name = os.path.join(folder, fname)
        (name_freq, name_to_id, prot_count) = \
                    fasta_lib.uniprot_species_frequency(db_name)
        #
        # sort the species names and write to file
        #
        fasta_lib.save_species_info(DB[i], folder, name_freq, \
                                    name_to_taxon, sci_to_taxon, \
                                    id_to_taxon, name_to_id, \
                                    minimum=min_sequence_count[i])
        #
        # print out some stats and exit
        #
        new_db = 'uniprot_%s_%s.fasta.gz' % (DB[i], versions[DB[i]],)
        for obj in write:
            print >>obj, '...%s contained %s entries...' % (new_db, prot_count)
            print >>obj, '...there were', len(name_freq), 'species names...'
    #
    fasta_lib.combine_analysis_files(folder)
    #
    fasta_lib.time_stamp_logfile('>>> ending: uniprot_get_analyze.py', log_obj)
    log_obj.close()
    #
    return
#
# end
#
if __name__ == '__main__':
    #
    import os
    import sys
    import fasta_lib
    #
    # get path to uniprot databases and call main function to download, etc.
    # check if folder path is passed on command line
    #
    versions = fasta_lib.get_uniprot_version()
    if len(sys.argv) > 1 and os.path.exists(sys.argv[1]):
        container = sys.argv[1]
    else:
        #
        # browse to a container folder for uniprot downloads
        # a subfolder will be created with name "uniprot_version"
        #
        database = r'C:\Xcalibur\database'  # default for BioWorks, XP
        if not os.path.exists(database):
            database = os.getcwd()
        container = fasta_lib.get_folder(database, \
                                         'Select folder for uniprot downloads')
        if container == '': sys.exit()    # cancel button response            
    folder = os.path.split(container)[1].lower()    # ignore case
    if folder.startswith('uniprot_'):   # if subfolder, up one level
        container = os.path.split(container)[0]
    folder = os.path.join(container, 'uniprot_' + versions['uniprot'])
    if not os.path.exists(folder):  # make folder if necessary
        os.mkdir(folder)
    #
    # pass in both databases for combined extraction
    #
    main(['sprot', 'trembl'], folder, versions)
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

                
    


