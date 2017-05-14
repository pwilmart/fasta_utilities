"""'check_for_duplicates.py' Written by Phil Wilmarth, OHSU.
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
    print '======================================================================='
    print ' check_for_duplicates.py, v1.01, written by Phil Wilmarth, OHSU, 2009.'
    print '======================================================================='
    #
    # set up the output file names, etc.
    #
    folder = os.path.split(fasta_file)[0]
    out_file = os.path.join(folder, 'duplicates.txt')
    out_obj = open(out_file, 'w')
    #
    # create instances of reader object and protein object, initialize counter
    #
    f = fasta_lib.FastaReader(fasta_file)
    p = fasta_lib.Protein()
    prot = 0
    head = 0
    dup = 0
    conflict = {}
    #
    # read proteins until EOF
    #
    while f.readNextProtein(p, check_for_errs=False):
        prot += 1
        control_A = p.description.count(chr(01))
        head = head + control_A + 1
        dup_data = (p.seqlenProtein(), p.molwtProtein())
        value = p.accession
        duplicate = conflict.get(dup_data, False)
        if duplicate:
            dup += 1
            print >>out_obj, '...WARNING: (protein no. %s) %s may be same as %s' % \
                  (prot, p.accession, duplicate)
        else:
            conflict[dup_data] = value
    #
    # print result and return
    #
    print '...there are %s proteins in %s' % \
          (prot, os.path.split(fasta_file)[1])
    if head > prot:
        print '...there were %s header lines...' % (head,)
    print '...there were %s possible duplicates...' % (dup,)
    out_obj.close()
    #
    # rewind out file and build new dictionaries
    #
    out_obj = open(out_file, 'r')
    conflict = {}
    to_save = {}
    while True:
        line = out_obj.readline()
        if not line:
            break
        else:
            line = line.split()
        conflict[line[4]] = line[9]
        to_save[line[4]] = True
        to_save[line[9]] = True
    out_obj.close()
    #
    # read in and save the proteins that might be duplicates of each other
    #
    candidates = []
    i = 0
    dup2 = 0
    index = {}
    f = fasta_lib.FastaReader(fasta_file)
    p = fasta_lib.Protein()
    while f.readNextProtein(p, check_for_errs=False):
        if to_save.get(p.accession, False):
            candidates.append(copy.deepcopy(p))
            index[p.accession] = i
            i += 1
    if len(candidates) == 0:
        print 'to_save_dictionary:', to_save
        print 'bail out in middle'
        return(candidates, index)
    #
    # look deeper to see if candidates are actually duplicates
    #
    out_obj = open(out_file, 'a')
    print >>out_obj, '\n========================================\n'
    exact_dup = 0
    accessions = conflict.keys()
    accessions.sort()
    for acc in accessions:
        p_dup = candidates[index[acc]]
        dup_acc = conflict[acc]
        p_ref = candidates[index[dup_acc]]
        if p_ref.sequence == p_dup.sequence:
            exact_dup += 1
            print >>out_obj, '...(%s) WARNING: %s is an exact match to %s' %\
                  (exact_dup, p_ref.accession, p_dup.accession)
    print '...number of exact matches was', exact_dup
    #
    return(candidates, index)
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
        database = r'E:\Ravi_thesis\Carr_plasma\databases\saved_DBs'
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
    candidates, index = main(fasta_file)
    try:    # wait for user to end program if not running from IDLE
        # what __file__ is under XP IDLE: 'C:\\Python26\\Lib\\idlelib\\idle.pyw'
        if not __file__.endswith('idle.pyw'):
            raw_input('\n...hit any key to end program...')
    except NameError:
        pass
#
# end
#
