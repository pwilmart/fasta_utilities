"""'nr_get_analyze.py' Written by Phil Wilmarth, OHSU.

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
# updated for Python 3 compatibility -PW 7/4/2017

import os
import sys
import time
import gzip
import fasta_lib

# minimum sequence count cutoff for output table
# NOTE: There are a LOT of species with less than 10 sequences
min_sequence_count = 10


def main(db, folder):
    """Fetches and analyzes the species names in the ncbi nr fasta database.

    Arguments:
    "nr_folder" is the full path name where "nr.gz" will be.  No return values.

    Saves the main lookup dictionary and a summary text file that
    can be loaded into EXCEL or a word processor.
    """
    global min_sequence_count

    print('================================================================')
    print(' nr_get_analyze.py, v1.1.0, written by Phil Wilmarth, OHSU, 2017')
    print('================================================================')

    # create a log file to mirror screen output
    log_obj = open(os.path.join(folder, 'fasta_utilities.log'), 'a')
    write = [None, log_obj]
    fasta_lib.time_stamp_logfile('\n>>> starting: nr_get_analyze.py', log_obj)

    # make sure the files are present or download if not
    fasta_lib.download_ncbi(folder)

    # make gi_to_taxon object (or reload from disk)
    acc_to_taxon = fasta_lib.AccToTaxon(folder)
    acc_to_taxon.create_or_load(folder)

    # make a dictionary of taxon IDs to species names
    taxon_to_name = fasta_lib.make_taxon_to_sci_name(folder)

    # make the taxon frequency dictionary for the proteins in nr.gz
    nr_name = os.path.split(folder)[1] + '.gz'
    for obj in write:
        print('...processing %s (this takes a few hours...)' % (nr_name,), file=obj)
    taxon_freq = {}
    reftax_freq = {}
    prot = 0
    spec_prot = 0
    ref_prot = 0
    undef_gi = 0
    fasta_file = gzip.open(os.path.join(folder, nr_name), mode='rt')
    while True:
        line = fasta_file.readline()
        if not line:
            break
        else:
            line = line.rstrip()
        if line.startswith('>'):
            prot += 1
            chunk = 1000
            if (prot % chunk) == 0:
                print('......(%s proteins read)' % ("{0:,d}".format(prot),))
            tax_list = []
            reftax_list = []
            for header in line[1:].split(chr(1)):   # need to remove the leading ">" character
                acc_ver = header.split()[0]
                acc = acc_ver.split('.')[0]
                tax = acc_to_taxon.get(acc, -1)
                if tax == -1:
                    undef_gi += 1
                if tax  not in tax_list:
                    spec_prot += 1
                    tax_list.append(tax)
                if '_' in acc and tax not in reftax_list:   # according to NCBI underscore char only in RefSeq
                    ref_prot +=1
                    reftax_list.append(tax)
            for tax in tax_list:
                fasta_lib.add_or_increment(tax, taxon_freq)
            for reftax in reftax_list:
                fasta_lib.add_or_increment(reftax, reftax_freq)

    # make the name frequency dictionary from the taxon frequency dictionary
    name_freq = {}
    for (taxon, freq) in taxon_freq.items():
        unknown_name = 'Unknown_taxonID_%s' % (taxon,)
        name_freq[taxon_to_name.get(taxon, unknown_name)] = freq

    # make an inverted name_to_taxon dictionary
    name_to_taxon = {}
    for (number, name) in taxon_to_name.items():
        name_to_taxon[name] = number

    # sort the species names and write to file
    fasta_lib.save_species_info_nr(folder, name_freq, name_to_taxon, reftax_freq, min_sequence_count)

    # print out some stats and exit
    for obj in write:
        if prot > 0:
            print('...%s contained %s protein entries...' %
                  (os.path.split(folder)[1], "{0:,d}".format(prot)), file=obj)
        print('...there were %s species-expanded entries...' % ("{0:,d}".format(spec_prot),), file=obj)
        print('...%s were RefSeq entries...' % ("{0:,d}".format(ref_prot),), file=obj)
        print('...%s entries had undefined taxon ID numbers...' % ("{0:,d}".format(undef_gi),), file=obj)
        print('...there were', "{0:,d}".format(len(name_freq)), 'species names...', file=obj)

    fasta_lib.time_stamp_logfile('>>> ending: nr_get_analyze.py', log_obj)
    log_obj.close()
    return


if __name__ == '__main__':
    # get the path to nr.gz and call main function to download, etc.

    # check if nr.gz file path is passed on command line
    if len(sys.argv) > 1 and os.path.exists(sys.argv[1]):
        selection = sys.argv[1]

    # if not, browse to a folder for nr downloads
    else:
        database = r'C:\Xcalibur\database'
        if not os.path.exists(database):
            database = os.getcwd()
        selection = fasta_lib.get_folder(database, 'Select folder for nr downloads')
        if selection == '': sys.exit() # cancel button response

    # if folder name starts with 'nr_', then skip creating a new folder
    if os.path.split(selection)[1].startswith('nr_'):
        main('nr', selection)

    # otherwise create a new folder with date stamp
    else:
        curr_time = time.localtime()
        curr_date = '%4d%02d%02d' % (curr_time[0], curr_time[1], curr_time[2])
        folder = os.path.join(selection, 'nr_'+curr_date)
        if not os.path.exists(folder):
            os.mkdir(folder)
        main('nr', folder)

# end
