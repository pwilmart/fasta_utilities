"""'count_fasta.py' Written by Phil Wilmarth, OHSU.
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
#=============================
def fasta_counter(fasta_file):
#=============================
    """Counts entries in a FASTA protein database.
        Call with FASTA filename, returns integer protein count
        Checks for duplicate accessions and (optional) valid characters.
    """
    import os
    import fasta_lib
    #
    print '==============================================================='
    print ' count_fasta.py, v 1.01, written by Phil Wilmarth, OHSU, 2009.'
    print '==============================================================='
    #
    # create a log file to mirror screen output
    #
    _folder = os.path.split(fasta_file)[0]
    log_obj = open(os.path.join(_folder, 'fasta_utilities.log'), 'a')
    write = [None, log_obj]
    fasta_lib.time_stamp_logfile('\n>>> starting: count_fasta.py', log_obj)
    #
    # create instances of reader object and protein object, initialize counters
    #
    f = fasta_lib.FastaReader(fasta_file)
    p = fasta_lib.Protein()
    prot = 0
    head = 0
    conflict = {}
    #
    # read proteins until EOF
    # NOTE: checking for errors slows program by factor of 3 or 4
    #
    while f.readNextProtein(p, check_for_errs=False):
        #
        # count protein sequences
        #
        prot += 1   
        if (prot % 500000) == 0:
            print '......(%s proteins read...)' % (prot,)
        #
        # check for duplicate accession
        #
        dup = conflict.get(p.accession, False)
        if dup:
            for obj in write:
                print >>obj, '\n...WARNING: %s is already in FASTA database!\n' % \
                      (p.accession,)
                if p.molwtProtein(show_errs=False) == conflict[p.accession]:
                    print >>obj, '......possible duplicated sequence...'
        else:
            conflict[p.accession] = p.molwtProtein(show_errs=False)
        #
        # count number of header elements
        #
        control_A = p.description.count(chr(01))
        head = head + control_A + 1
    #
    # print results and return
    #
    for obj in write:
        print >>obj, '...there are %s proteins in %s' % \
              (prot, os.path.split(fasta_file)[1])
        if head > prot:
            print >>obj, '...there were %s header lines' % (head,)
    #
    fasta_lib.time_stamp_logfile('>>> ending: count_fasta.py', log_obj)
    log_obj.close()
    #    
    return(prot)
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
    # call counter function
    #
    fasta_counter(fasta_file)
    try:    # wait for user to end program if not running from IDLE
        # what __file__ is under XP IDLE: 'C:\\Python26\\Lib\\idlelib\\idle.pyw'
        if not __file__.endswith('idle.pyw'):
            raw_input('\n...hit any key to end program...')
    except NameError:
        pass
#
# end
#
