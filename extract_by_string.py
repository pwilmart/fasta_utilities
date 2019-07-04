"""'extract_by_string.py' Written by Phil Wilmarth, OHSU.

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
import fasta_lib

# flag for case-sensitive matching (True) or not (False)
CASE_SENSITIVE = True

# clean accessions/descriptions or not (and DB-specific cleaning options)
# NOTE: information can be lost if accessions/descriptions are cleaned
CLEAN_ACCESSIONS = False
REF_SEQ_ONLY = True
KEEP_UNIPROT_ID = False

# list strings to search for in header lines and name to use in filenames (do not include leading underscores)
# left string is the search string, right string is the output file tag (underscore added sutomatically)
##string_dict = { '|ref|':'ref_seq' }
##string_dict = { '_MYCLB':'_MYCLB', '_MYCLE':'_MYCLE' }
##string_dict = { '[Bacillus sp.':'Bacillus_sp' }
string_dict = { 'Homo sapiens':'human_string'}
string_dict = { 'Uncharacterized protein':'uncharacterized_proteins'}


def main(string_dict):
    """Main program to extract entries containing strings from databases.
        Simple string search of pattern in combined accession/description lines.
        Logical OR if more than one pattern is mapped to the same outfile.
        Each matching protein is written once per output file with possible
            compound header (nr) of all headers containing matching patterns.
            If "cleaning" of accessions/descriptions is turned on for NCBI nr
            databases, only the first header element will be retained and
            any accession number cross-references will be lost.

    Written by Phil Wilmarth, OHSU, 2009.
    """
    print('=====================================================================')
    print(' extract_by_string.py, v.1.1.0, written by Phil Wilmarth, OHSU, 2017 ')
    print('=====================================================================')

    # set some file paths and names
    default = r'C:\Xcalibur\database'
    if not os.path.exists(default):
        default = os.getcwd()
    db_file = fasta_lib.get_file(default, [('Zipped files', '*.gz'), ('Fasta files', '*.fasta')],
                                 title_string='Select a FASTSA database')
    if db_file == '' : sys.exit() # cancel button repsonse

    db_folder, db_name = os.path.split(db_file)
    base_name = db_name.replace('.gz', '')
    if not base_name.endswith('.fasta'):
        base_name = base_name + '.fasta'

    # create a log file to mirror screen output
    log_obj = open(os.path.join(db_folder, 'fasta_utilities.log'), 'a')
    write = [None, log_obj]
    fasta_lib.time_stamp_logfile('\n>>> starting: extract_by_string.py', log_obj)

    # print the list of patterns that will be extracted
    string_list = list(string_dict.items())
    string_list.sort()
    for obj in write:
        print('...extracting entries containing these strings:', file=obj)
        for i, t in enumerate(string_list):
            print('......(%s) string "%s" to file ending in "%s"' % (i+1, t[0], t[1]), file=obj)

    # open the output databases, initialize counters
    string_files = {}
    string_count = {}
    name_count = {}
    for string, name in string_dict.items():
        fname = base_name.replace('.fasta', '_'+name+'.fasta')
        fname = os.path.join(db_folder, fname)
        string_files[name] = fname
        string_count[string] = 0
        name_count[name] = 0
    for name in string_files.keys():
        string_files[name] = open(string_files[name], 'w')

    # create a FastaReader object, initialize counters, and start reading
    x = fasta_lib.FastaReader(db_file)
    prot = fasta_lib.Protein()
    prot_read = 0
    for obj in write:
        print('...reading %s and extracting entries...' % (db_name,), file=obj)
    while x.readNextProtein(prot, check_for_errs=False):
        prot_read += 1
        if (prot_read % 500000) == 0:
            print('......(%s proteins read...)' % ("{0:,d}".format(prot_read),))
        written = {}    # make sure protein is written only ONCE per OUTFILE
        header = prot.accession + ' ' + prot.description # recreate the '>' line
        if not CASE_SENSITIVE:  # convert to uppercase
            header = header.upper()
        for pattern in string_dict.keys():
            new_pattern = pattern
            if not CASE_SENSITIVE:  # case insensitive matching
                new_pattern = new_pattern.upper()
            for head in header.split(chr(1)):  # check each header for matches
                if new_pattern in head:
                    name = string_dict[pattern]
                    name_header = written.get(name, '')
                    if name_header:
                        name_header = name_header + chr(1) + head
                        written[name] = name_header
                    else:
                        written[name] = head
                        string_count[pattern] += 1

        # write any matching proteins to appropriate out file
        for name in written.keys():
            name_count[name] += 1   # output file write counters
            f = string_files[name]  # output file pointers
            header = written[name]  # composite header of name's matches

            # set the accession and description fields before writing
            prot.accession = header.split()[0]
            prot.new_acc = prot.accession
            prot.description = header[(len(prot.accession)+1):]
            prot.new_desc = prot.description
            if CLEAN_ACCESSIONS:
                if prot.accession.startswith('gi|'):
                    prot.parseNCBI(REF_SEQ_ONLY)
                elif prot.accession.startswith('sp|') or prot.accession.startswith('tr|'):
                    prot.parseUniProt(KEEP_UNIPROT_ID)
            prot.printProtein(f)    # write any matching proteins

    # close files
    for f in string_files.values():
        f.close()

    # print out the summary stuff
    for obj in write:
        print('...%s protein entries in %s' % ("{0:,d}".format(prot_read), db_name), file=obj)
        strings = list(string_count.keys())
        strings.sort()
        for i, string in enumerate(strings):
            print('......(%s) pattern "%s" was found in %s proteins' %
                  (i+1, string, "{0:,d}".format(string_count[string])), file=obj)
        print('...output file summaries...', file=obj)
        names = list(string_files.keys())
        names.sort()
        for i, name in enumerate(names):
            temp = base_name.replace('.fasta', '_'+name+'.fasta')
            print('......(%s) %s proteins extracted and written to %s' %
                  (i+1, "{0:,d}".format(name_count[name]), temp), file=obj)

    fasta_lib.time_stamp_logfile('>>> ending: extract_by_string.py', log_obj)
    log_obj.close()
    return

# check for command line launch and see if any arguments passed
if __name__ == '__main__':
    if len(sys.argv) > 1:
        arg_dict = fasta_lib.string_cmd_line_checker(sys.argv)
        if arg_dict:
            main(arg_dict)
        else:
            sys.exit()
    else:
        main(string_dict)

# end
