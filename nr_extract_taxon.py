"""'nr_extract_taxon.py' Written by Phil Wilmarth, OHSU.
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
EXPAND_GROUPS = False        # controls expanding taxonomy groups (nodes)
MIN_GROUP_SEQ_COUNT = 10     # minimum if expanding a taxon group
#
# extract only RefSeq entries if True
#
REF_SEQ_ONLY = True
#
# other flags (NOTE: information can be lost by cleaning accessions)
#
CLEAN_ACCESSIONS = False    # one header element, either gi or RefSeq accession
VERBOSE = True              # prints more information for taxon nodes
#
# list species to extract by taxonomy number and name to use in filenames
#
taxon_dict = { 9606:'human_refseq',
               10090:'mouse_refseq',
               10116:'rat_refseq',
               559292:'yeast_refseq',
               9544:'Macaca.mulatta',
               83332: 'M_tuberculosis_H37Rv' } # default list of species
##taxon_dict = { 9606:'human_all',
##               10090:'mouse_all',
##               10116:'rat_all',
##               559292:'yeast_all'} # non-RefSeq
##taxon_dict = { 7108: 'S_frugiperda',
##               7091: 'B_mori' } # insects
#
#
#====================
def main(taxon_dict):
#====================
    """Main program to extract entries by taxon ID from NCBI nr databases.
        Each gi number (of each header) is looked up to find associated taxon
        number for comparison to desired taxon numbers.  A separate protein
        entry will be written for each desired taxon number even if all taxon
        numbers are written to the same output file.  At the protein level, the
        extracted databases may no longer be non-redundant.  If "cleaning" of
        accessions/descriptions is turned off, all headers matching the desired
        taxon numbers will be added to the respective protein preserving the
        usual NCBI nr formatting structure.  If cleaning of accessions is turned
        on during extraction, some information may be lost.  This could make
        subsequent database processing (such as extracting by text string) fail.
        Cleaning is best done as a last step (i.e. in "reverse_fasta.py").
    """
    import fasta_lib
    import cPickle
    import os
    import sys
    import copy
    #
    print '===================================================================='
    print ' nr_extract_taxon.py, v.1.02, written by Phil Wilmarth, OHSU, 2009.'
    print '===================================================================='
    #
    # set some file paths and names
    #
    default = r'C:\Xcalibur\database'
    if not os.path.exists(default):
        default = os.getcwd()
    nr_file = fasta_lib.get_file(default, [('Zipped files', '*.gz'), \
                                           ('Fasta files', '*.fasta')], \
                                 title_string='Select an NCBI nr database')
    if nr_file == '': sys.exit() # cancel button response
    ncbi_folder, nr_name = os.path.split(nr_file)
    nr_db = os.path.splitext(nr_name)[0]
    #
    # create a log file to mirror screen output
    #
    log_obj = open(os.path.join(ncbi_folder, 'fasta_utilities.log'), 'a')
    write = [None, log_obj]
    fasta_lib.time_stamp_logfile('\n>>> starting: nr_extract_taxon.py', log_obj)
    #
    # get the saved gi number to taxon number {int:int} dictionary
    #
    gi_to_taxon = fasta_lib.GiToTaxon(ncbi_folder)
    gi_to_taxon.create_or_load(ncbi_folder)
    #
    # print the list of taxon numbers that will be extracted
    #
    original_dict = taxon_dict
    taxon_list = taxon_dict.items()
    taxon_list.sort()
    for obj in write:
        print >>obj, '...extracting these taxon numbers:'
        for i, t in enumerate(taxon_list):
            print >>obj, '......(%s) taxon %s to file tagged with "%s"' % \
                  (i+1, t[0], t[1])
    #
    # expand any group taxon numbers.  NOTE: if a taxon number appears in
    # "nr_fasta_analyze.txt", it will not be expanded.  Either delete the
    # line in "nr_fasta_analyze.txt", or make an expanded "taxon_dict" by hand.
    #
    if EXPAND_GROUPS:
        fasta_lib.expand_species(ncbi_folder, 'nr', taxon_dict, \
                                 MIN_SEQUENCE_COUNT, \
                                 MIN_GROUP_SEQ_COUNT, REF_SEQ_ONLY)
    #
    # open the output databases, initialize counters, etc.
    #
    taxon_files = {}
    taxon_count = {}
    name_count = {}
    for taxon, name in taxon_dict.items():
        fname = nr_db+'_'+name+'.fasta'
        fname = os.path.join(ncbi_folder, fname)
        taxon_files[name] = fname
        name_count[name] = 0
        taxon_count[taxon] = 0
    #
    # open the output filenames
    #
    for name in taxon_files.keys():
        taxon_files[name] = open(taxon_files[name], 'w')
    #
    # loop over all proteins in nr
    #
    x = fasta_lib.FastaReader(nr_file)
    prot = fasta_lib.Protein()
    prot_read = 0
    not_found = 0
    skipped = 0
    for obj in write:
        print >>obj, '...reading %s and extracting entries...' % (nr_name,)
    #
    # checking for errors slows down program by about a factor of 3 or 4
    #
    while x.readNextProtein(prot, check_for_errs=False):
        prot_read += 1
        if (prot_read % 500000) == 0:
            print '......(%s proteins read...)' % (prot_read,)
        written = {}
        line = prot.accession + ' ' + prot.description
        prot.new_desc = ''
        #
        # extract the gi numbers for each header
        #
        for header in line.split(chr(01)):
            accession = header.split()[0]
            if REF_SEQ_ONLY and '|ref|' not in accession:
                continue    # skip proteins without RefSeq entries
            gi_number = int(accession.split('|')[1])
            taxon = gi_to_taxon.get(gi_number, False)
            #
            # see if taxon number for this gi is in our desired list
            #
            if taxon:
                if taxon_dict.get(taxon, False):
                    if written.get(taxon, False):
                        # if taxon number already seen, add to header
                        prot = written[taxon]
                        prot.description = prot.description + chr(01) + header
                        written[taxon] = copy.deepcopy(prot)
                    else:
                        # first time taxon number seen
                        name = taxon_dict[taxon]
                        prot.accession = header.split()[0]
                        prot.description = header[len(prot.accession)+1:]
                        prot.description = prot.description.rstrip()
                        taxon_count[taxon] += 1
                        name_count[name] += 1
                        written[taxon] = copy.deepcopy(prot)
                else:
                    skipped += 1
            else:
                not_found += 1
                continue
        #
        # write a protein sequence for each taxon number it was matched to
        #        
        for taxon in written.keys():
            name = taxon_dict[taxon]
            f = taxon_files[name]
            prot = written[taxon]
            prot.new_desc = prot.description
            prot.new_acc = prot.accession
            if CLEAN_ACCESSIONS:
                prot.parseNCBI(REF_SEQ_ONLY)
            prot.printProtein(f)
    #
    # print out number of matches and close files
    #
    for obj in write:
        print >>obj, '...%s proteins in %s' % (prot_read, nr_name)
        print >>obj, '...%s gi numbers did not have known taxon numbers' % \
              (not_found,)
        print >>obj, '...%s gi numbers did not have matching taxon numbers' % \
              (skipped,)
        if REF_SEQ_ONLY:
            print >>obj, '...Extracted sequences are RefSeq Only!!!'
        if VERBOSE:
            numbers = taxon_count.keys()
            numbers.sort()
            for i, number in enumerate(numbers):
                if taxon_count[number] > 0:
                    print >>obj, '......(%s) taxon number %s had %s proteins' % \
                          (i+1, number, taxon_count[number])
        print >>obj, '...output file summaries...'
        names = taxon_files.keys()
        names.sort()
        for i, name in enumerate(names):
            print >>obj, '......(%s) %s proteins extracted and written to %s' % \
                  (i+1, name_count[name], nr_db+'_'+name+'.fasta')
    #
    fasta_lib.time_stamp_logfile('>>> ending: nr_extract_taxon.py', log_obj)
    log_obj.close()
    for f in taxon_files.values():
        f.close()
    #
    return
#
# check for command line launch and see if any arguments passed
#
import sys
import os
import fasta_lib
if __name__ == '__main__':
    #
    # if arguments make sure they are taxon name pairs
    #
    if len(sys.argv) > 1:
        arg_dict = fasta_lib.taxon_cmd_line_checker(sys.argv)
        if arg_dict:
            main(arg_dict)
        else:
            sys.exit()  # error in command line arguments
    else:
        main(taxon_dict)
    try:    # wait for user to end program if not running from IDLE
        # what __file__ is under XP IDLE: 'C:\\Python26\\Lib\\idlelib\\idle.pyw'
        if not __file__.endswith('idle.pyw'):
            raw_input('\n...hit any key to end program...')
    except NameError:
        pass
#
# end
#
    



         
