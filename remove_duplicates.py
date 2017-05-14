"""'remove_duplicates.py' Written by Phil Wilmarth, OHSU.
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
#
#================================================================
def find_identities(prot, candidates, skip, fasta_file, out_obj):
#================================================================
    """Test sequences in candiates list for identity to prot.
    """
    ident = 0
    acc = prot.accession
    seq = prot.sequence
    for p in candidates[acc]:
        if seq == p.sequence:
            ident += 1
            skip[p.accession] = True
            prot.new_desc = prot.new_desc + chr(01) + p.accession + ' ' + p.description
            if ident == 1:
                print >>out_obj, '\n...%s "%s" same as:' % (acc, prot.description[:60])
                print >>out_obj, '......%s "%s"' % (p.accession, p.description[:60])
            else:
                print >>out_obj, '......%s "%s"' % (p.accession, p.description[:60])
    return ident
#
#
#====================
def main(fasta_file):
#====================
    """Checks entries in a FASTA protein database for identical duplicates.
        Call with FASTA filename, returns a couple of dictionaries
    """
    import os
    import fasta_lib
    import copy
    #
    # set up the output file names, etc.
    #
    folder = os.path.split(fasta_file)[0]
    out_file = os.path.join(folder, 'duplicates.txt')
    out_obj = open(out_file, 'a')
    for end in ['.fasta', '.fasta.gz', '.fa.gz']:
        if fasta_file.endswith(end):
            nr_database = fasta_file.replace(end, '_nonredun.fasta')
    if (not nr_database) or (nr_database == fasta_file):
        nr_database = fasta_file + '_nonredun.fasta'
    nr_obj = open(nr_database, 'w')
    write = [None, out_obj]
    fasta_lib.time_stamp_logfile('\n>>> starting: check_for_duplicates.py', out_obj)
    #
    print '===================================================================='
    print ' remove_duplicates.py, v1.02, written by Phil Wilmarth, OHSU, 2009. '
    print '===================================================================='
    print '\n   Analysis will be written to "duplicates.txt"'
    #
    # create instances of reader object and protein object, initialize counter
    #
    f = fasta_lib.FastaReader(fasta_file)
    p = fasta_lib.Protein()
    prot = 0
    head = 0
    dup = 0
    candidates = {}     # dictionary of acc:protein lists
    conflicts = {}      # keeps track of seq len and MW
    #
    # read proteins until EOF
    #
    while f.readNextProtein(p, check_for_errs=False):
        prot += 1
        control_A = p.description.count(chr(01))
        head = head + control_A + 1
        dup_data = (p.seqlenProtein(), p.molwtProtein())
        duplicate = conflicts.get(dup_data, False)
        if duplicate:
            dup += 1
            if candidates.get(duplicate, False):
                candidates[duplicate].append(copy.deepcopy(p))
            else:
                candidates[duplicate] = [copy.deepcopy(p)]
        else:
            conflicts[dup_data] = p.accession
    #
    # get list of proteins to test for identity
    #
    to_test = {}
    for key in candidates.keys():
        to_test[key] = True    
    #
    # copy proteins to "nr" file, checking for duplicates
    #
    dup = 0
    skip = {}
    f = fasta_lib.FastaReader(fasta_file)
    p = fasta_lib.Protein()
    print >>out_obj, 'Processing:', fasta_file    # header line to log file
    while f.readNextProtein(p, check_for_errs=False):
        if to_test.get(p.accession, False):
            dup += find_identities(p, candidates, skip, fasta_file, out_obj)
        if skip.get(p.accession, False):
            continue
        p.printProtein(nr_obj)
    #
    for obj in [None, out_obj]:
        print >>obj, '\nThere were', prot, 'total sequences in:', os.path.basename(fasta_file)
        print >>obj, 'There were', dup, 'identical sequences removed\n\n'
        try:
            obj.close()
        except AttributeError:
            pass
    return
    #
    # end
    #
#
# setup stuff: check for command line args, etc.
#
if __name__ == '__main__':
    import os
    import sys
    import fasta_lib
    #
    # check if database name passed on command line
    #
    if len(sys.argv) > 1 and os.path.exists(sys.argv[1]):
        fasta_file = sys.argv[1]
    #
    # if not, browse to database file
    #    
    else:
        if len(sys.argv) > 1:
            print '...WARNING: %s not found...' % (sys.argv[1],)
        database = r'C:\Xcalibur\database'
        if not os.path.exists(database):
            database = os.getcwd()
        fasta_file = fasta_lib.get_file(database, \
                                        [('FASTA files', '*.fasta'),\
                                         ('Zipped FASTA files', '*.gz'),\
                                         ('All files', '*.*')],\
                                        'Select a FASTA database')
        if fasta_file == '': sys.exit()     # cancel button repsonse
    #
    # call main function
    #
    main(fasta_file)
    try:    # wait for user to end program if not running from IDLE
        # what __file__ is under XP IDLE: 'C:\\Python26\\Lib\\idlelib\\idle.pyw'
        if not __file__.endswith('idle.pyw'):
            raw_input('\n...hit any key to end program...')
    except NameError:
        pass
#
# end
#
