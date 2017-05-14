"""'ipi_parser.py' Written by Phil Wilmarth, OHSU.
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
# flag to keep the Gene identifier in the description
#
KEEP_IPI_GENE_ID = True
#
#==========================
def ipi_parser(fasta_file):
#==========================
    """Parses accessions and descriptions in IPI databases.
        Call with FASTA filename, no return values.
    """
    import os
    import fasta_lib
    #
    print '=============================================================='
    print ' ipi_parser.py, v 1.01, written by Phil Wilmarth, OHSU, 2009.'
    print '=============================================================='
    #
    # create a log file to mirror screen output
    #
    _folder = os.path.split(fasta_file)[0]
    log_obj = open(os.path.join(_folder, 'fasta_utilities.log'), 'a')
    write = [None, log_obj]
    fasta_lib.time_stamp_logfile('\n>>> starting: ipi_parser.py', log_obj)
    #
    # create instances of reader object and protein object, initialize counter
    #
    out_name = fasta_file.replace('.fasta', '.clean.fasta')
    if out_name.endswith('.gz'):
        out_name = out_name[:-3]
    out_file_obj = open(out_name, 'w')
    f = fasta_lib.FastaReader(fasta_file)
    p = fasta_lib.Protein()
    prot = 0
    notipi = 0
    #
    # read proteins and parse accessions, descriptions until EOF
    # NOTE: error checking slows program execution by factor of 3 or 4
    #
    while f.readNextProtein(p, check_for_errs=False):
        prot += 1
        if p.accession.startswith('IPI:'):
            p.parseIPI(KEEP_IPI_GENE_ID)
            p.printProtein(out_file_obj)
        else:
            notipi += 1
    #
    # print result and return
    #
    for obj in write:
        print >>obj, '...there were %s proteins cleaned in %s' % \
              (prot, os.path.split(fasta_file)[1])
        if notipi:
            print >>obj, '\n...WARNING: %s proteins did not have IPI formatting' \
                  % (notipi,)
    #
    fasta_lib.time_stamp_logfile('>>> ending: ipi_parser.py', log_obj)
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
                                        [('Zipped FASTA files', '*.gz'),\
                                         ('FASTA files', '*.fasta'),\
                                         ('All files', '*.*')],\
                                        'Select a FASTA database')
        if fasta_file == '': sys.exit() # cancel button response
    #
    # call ipi_parser function
    #
    ipi_parser(fasta_file)
    try:    # wait for user to end program if not running from IDLE
        # what __file__ is under XP IDLE: 'C:\\Python26\\Lib\\idlelib\\idle.pyw'
        if not __file__.endswith('idle.pyw'):
            raw_input('\n...hit any key to end program...')
    except NameError:
        pass
#
# end
#
