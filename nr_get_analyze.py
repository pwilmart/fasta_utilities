"""'nr_get_analyze.py' Written by Phil Wilmarth, OHSU.
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
# minimum sequence count cutoff for output table
# NOTE: There are a LOT of species with less than 10 sequences
#
min_sequence_count = 10
#
#
#====================
def main(db, folder):
#====================
    """Fetches and analyzes the species names in the ncbi nr fasta database.

    Arguments:
    "nr_folder" is the full path name where "nr.gz" will be.  No return values.

    Saves the main lookup dictionary and a summary text file that
    can be loaded into EXCEL or a word processor.
    """
    import os
    import gzip
    import cPickle
    import fasta_lib
    global min_sequence_count
    #
    print '================================================================='
    print ' nr_get_analyze.py, v1.01, written by Phil Wilmarth, OHSU, 2009.'
    print '================================================================='
    #
    # create a log file to mirror screen output
    #
    log_obj = open(os.path.join(folder, 'fasta_utilities.log'), 'a')
    write = [None, log_obj]
    fasta_lib.time_stamp_logfile('\n>>> starting: nr_get_analyze.py', log_obj)
    #
    # make sure the files are present or download if not
    #
    fasta_lib.download_ncbi(folder)
    #
    # make gi_to_taxon object (or reload from disk)
    #
    gi_to_taxon = fasta_lib.GiToTaxon(folder)
    gi_to_taxon.create_or_load(folder)
    #
    # make a dictionary of taxon IDs to species names
    #
    taxon_to_name = fasta_lib.make_taxon_to_sci_name(folder)
    #
    # make the taxon frequency dictionary for the proteins in nr.gz
    #
    nr_name = os.path.split(folder)[1] + '.gz'
    for obj in write:
        print >>obj, '...processing %s (this can take many minutes...)' % \
              (nr_name,)
    taxon_freq = {}
    reftax_freq = {}
    prot = 0
    spec_prot = 0
    ref_prot = 0
    undef_gi = 0
    fasta_file = gzip.open(os.path.join(folder, nr_name), 'rb')
    while True:
        line = fasta_file.readline()
        if not line:
            break
        else:
            line = line.rstrip()
        if line.startswith('>'):
            prot += 1
            if (prot % 500000) == 0:
                print '......(%s proteins read)' % (prot,)
            tax_list = []
            reftax_list = []
            for header in line.split(chr(01)):
                gi = header.split()[0]
                gi = int(gi.split('|')[1])
                tax = gi_to_taxon.get(gi, -1)
                if tax == -1:
                    undef_gi += 1
                if tax  not in tax_list:
                    spec_prot += 1
                    tax_list.append(tax)
                if '|ref|' in header and tax not in reftax_list:
                    ref_prot +=1
                    reftax_list.append(tax)
            for tax in tax_list:
                fasta_lib.add_or_increment(tax, taxon_freq)
            for reftax in reftax_list:
                fasta_lib.add_or_increment(reftax, reftax_freq)                
    #
    # make the name frequency dictionary from the taxon frequency dictionary
    #
    name_freq = {}
    for (taxon, freq) in taxon_freq.items():
        unknown_name = 'Unknown_taxonID_%s' % (taxon,)
        name_freq[taxon_to_name.get(taxon, unknown_name)] = freq
    #
    # make an inverted name_to_taxon dictionary
    #
    name_to_taxon = {}
    for (number, name) in taxon_to_name.items():
        name_to_taxon[name] = number
    #
    # sort the species names and write to file
    #
    fasta_lib.save_species_info_nr(folder, name_freq, name_to_taxon, \
                                     ref2freq=reftax_freq, \
                                     minimum=min_sequence_count)
    #
    # print out some stats and exit
    #
    for obj in write:
        if prot > 0:
            print >>obj, '...%s contained %s protein entries...' %\
                  (os.path.split(folder)[1], prot)
        print >>obj, '...there were %s species-expanded entries...' %\
              (spec_prot,)
        print >>obj, '...%s were RefSeq entries...' % (ref_prot,)
        print >>obj, '...%s entries had undefined taxon ID numbers...' % \
              (undef_gi,)
        print >>obj, '...there were', len(name_freq), 'species names...'
    #
    fasta_lib.time_stamp_logfile('>>> ending: nr_get_analyze.py', log_obj)
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
    import time
    import fasta_lib
    #
    # get the path to nr.gz and call main function to download, etc.
    # check if nr.gz file path is passed on command line
    #
    if len(sys.argv) > 1 and os.path.exists(sys.argv[1]):
        selection = sys.argv[1]
    else:
        #
        # browse to a folder for nr downloads
        #
        database = r'C:\Xcalibur\database'
        if not os.path.exists(database):
            database = os.getcwd()
        selection = fasta_lib.get_folder(database, \
                                         'Select folder for nr downloads')
        if selection == '': sys.exit() # cancel button response
    #
    # if folder name starts with 'nr_', then skip creating a new folder
    #
    if os.path.split(selection)[1].startswith('nr_'):
        main('nr', selection)
    #
    # otherwise create a new folder with date stamp
    #
    else:
        curr_time = time.localtime()
        curr_date = '%4d%02d%02d' % (curr_time[0], curr_time[1], curr_time[2])
        folder = os.path.join(selection, 'nr_'+curr_date)
        if not os.path.exists(folder):
            os.mkdir(folder)
        main('nr', folder)
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


                
    


