"""'uniprot_get_analyze.py' Written by Phil Wilmarth, OHSU.

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

# minimum sequence count cutoff for output table -
# there are MANY species in Trembl with 5 or fewer protein sequences
min_sequence_count = [0, 10]     # [sprot, trembl]


def main(DB, folder, versions):
    """Analyzes the species names in both UniProt databases.

    Arguments:
    "folder" is full path name to UniProt DBs.  "versions" is
    a dictionary of version numbers for file naming.  No return values.

    Saves a summary text file that can be loaded into EXCEL.
    """
    print('======================================================================')
    print(' uniprot_get_analyze.py, v1.1.0, written by Phil Wilmarth, OHSU, 2017 ')
    print('======================================================================')
    
    # create a log file to mirror screen output
    log_obj = open(os.path.join(folder, 'fasta_utilities.log'), 'a')
    write = [None, log_obj]
    fasta_lib.time_stamp_logfile('\n>>> starting: uniprot_get_analyze.py', log_obj)
    
    # make sure the files are present or download if not
    for db in DB:
        fasta_lib.download_uniprot(db, folder, versions)
    
    # make dictionaries of species names (or accession IDs) to taxon IDs
    (sci_to_taxon, id_to_taxon) = fasta_lib.make_uniprot_to_taxon(folder)
    
    # get more complete list of names to taxon number from ncbi data
    name_to_taxon = fasta_lib.make_all_names_to_taxon(folder)
    
    # make species frequency dictionary
    for i in range(len(DB)):
        fname = 'uniprot_%s_%s.fasta.gz' % (DB[i], versions[DB[i]],)
        db_name = os.path.join(folder, fname)
        (name_freq, name_to_id, prot_count) = fasta_lib.uniprot_species_frequency(db_name)
        
        # sort the species names and write to file
        fasta_lib.save_species_info(DB[i], folder, name_freq, name_to_taxon, sci_to_taxon,
                                    id_to_taxon, name_to_id, minimum=min_sequence_count[i])
        
        # print out some stats and exit
        new_db = 'uniprot_%s_%s.fasta.gz' % (DB[i], versions[DB[i]],)
        for obj in write:
            print('...%s contained %s entries...' % (new_db, prot_count), file=obj)
            print('...there were', len(name_freq), 'species names...', file=obj)
    
    fasta_lib.combine_analysis_files(folder)
    
    fasta_lib.time_stamp_logfile('>>> ending: uniprot_get_analyze.py', log_obj)
    log_obj.close()
    return


if __name__ == '__main__':
    # get path to uniprot databases and call main function to download, etc.
    # check if folder path is passed on command line
    versions = fasta_lib.get_uniprot_version()
    if len(sys.argv) > 1 and os.path.exists(sys.argv[1]):
        container = sys.argv[1]
    else:
        
        # browse to a container folder for uniprot downloads
        # a subfolder will be created with name "uniprot_version"
        database = r'C:\Xcalibur\database'  # default for BioWorks, XP
        if not os.path.exists(database):
            database = os.getcwd()
        container = fasta_lib.get_folder(database, 'Select folder for uniprot downloads')
        if container == '': sys.exit()    # cancel button response
        
    folder = os.path.split(container)[1].lower()    # ignore case
    if folder.startswith('uniprot_'):   # if subfolder, up one level
        container = os.path.split(container)[0]
    folder = os.path.join(container, 'uniprot_' + versions['uniprot'])
    if not os.path.exists(folder):  # make folder if necessary
        os.mkdir(folder)
    
    # pass in both databases for combined extraction
    main(['sprot', 'trembl'], folder, versions)

# end


                
    


