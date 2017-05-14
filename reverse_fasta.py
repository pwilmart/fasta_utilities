"""'reverse_fasta.py' Written by Phil Wilmarth, OHSU.
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
# flag to parse (clean) accessions and descriptions.
# NOTE: cleaning can cause some loss of information.  Apply at later stages
# of database processing such as here.
#
CLEAN_ACCESSIONS = False
#
# flag to keep UniProt Identifier (True) or more-stable ACC (False)
#
KEEP_UNIPROT_ID = False
#
# flag set if RefSeq entries are being extracted from NCBI nr
#
REF_SEQ_ONLY = True
#
# flag to keep IPI gene identifiers (True) or not (False)
#
KEEP_IPI_GENE_ID = True
#
# flags to make different output databases
#
MAKE_SEPARATE_FORWARD = True
MAKE_SEPARATE_REVERSED = False
MAKE_SEPARATE_BOTH = True
#
#=============================
def fasta_reverse(fasta_file):
#=============================
    """Adds contaminants and reverses entries in a FASTA protein database.

    Called with FASTA filename.  Reversed DB written to same location.
    Options for separate or concatenated output files.
    """
    import os
    decoy_string = 'REV_'   # the string to denote decoy sequences
    #
    print '================================================================'
    print ' reverse_fasta.py, v1.01, written by Phil Wilmarth, OHSU, 2009.'
    print '================================================================'
    #
    # open the "forward" and "reversed" output files
    #
    if fasta_file.lower().endswith('.gz'):
        _file = os.path.splitext(fasta_file[:-3])[0]
    else:
        _file = os.path.splitext(fasta_file)[0]
    for_name = _file + '_for.fasta'
    for_file_obj = open(for_name, 'w')
    rev_name = _file + '_rev.fasta'
    rev_file_obj = open(rev_name, 'w')
    #
    # create a log file to mirror screen output
    #
    _folder = os.path.split(fasta_file)[0]
    log_obj = open(os.path.join(_folder, 'fasta_utilities.log'), 'a')
    write = [None, log_obj]
    fasta_lib.time_stamp_logfile('\n>>> starting: reverse_fasta.py', log_obj)
    #
    # create instances of reader object and protein object
    # add the contaminants first
    #
    prot = fasta_lib.Protein()
    pcount = 0
    print 'cwd:', os.getcwd()
    try:
        if os.path.exists('all_contams_fixed.fasta'):
            print 'contams file found:', _file
            _file = 'all_contams_fixed.fasta'
        else:
            path = os.path.split(fasta_file)[0]
            _file = os.path.join(path, 'all_contams_fixed.fasta')
            print 'trying:', _file
        f = fasta_lib.FastaReader(_file)
        while f.readNextProtein(prot, check_for_errs=True):
            pcount += 1
            if CLEAN_ACCESSIONS:
                prot.parseCONT()
            prot.printProtein(for_file_obj)
            rev = prot.reverseProtein(decoy_string)
            rev.printProtein(rev_file_obj)
        for obj in write:
            print >>obj, '...there were %s contaminant entries in %s' % \
                  (pcount, _file)            
    except:
        for obj in write:
            print >>obj, '   WARNING: "all_contams_fixed.fasta" not found!'
    #
    # read proteins, clean up accessions, decriptions until EOF
    # write proteins to "forward" and "reversed" files
    #
    f = fasta_lib.FastaReader(fasta_file)
    #
    # error checking slows program execution, turn on if needed.
    # Reading and writing sequences always removes spaces and blank lines.
    #
    while f.readNextProtein(prot, check_for_errs=False):
        pcount += 1
        if CLEAN_ACCESSIONS:
            if prot.accession.startswith('gi|'):
                prot.parseNCBI(REF_SEQ_ONLY)
            elif prot.accession.startswith('sp|') or \
                 prot.accession.startswith('tr|'):
                prot.parseUniProt(KEEP_UNIPROT_ID)
            elif prot.accession.startswith('IPI:'):
                prot.parseIPI(KEEP_IPI_GENE_ID)
            else:
                pass
        #
        prot.printProtein(for_file_obj)    # write to "forward" file
        rev = prot.reverseProtein(decoy_string)
        rev.printProtein(rev_file_obj)   # write to "reversed" file
    #
    # make concatenated output file if desired and print summary stats
    #
    if MAKE_SEPARATE_BOTH:
        _file = fasta_file.replace('.gz', '')
        both_name = _file.replace('.fasta', '_both.fasta')
        both_file_obj = open(both_name, 'w')
        for_file_obj.close()
        for_file_obj = open(for_name, 'r')
        rev_file_obj.close()
        rev_file_obj = open(rev_name, 'r')
        while True:
            line = for_file_obj.readline()
            if not line: break
            both_file_obj.write(str(line))
        while True:
            line = rev_file_obj.readline()
            if not line: break
            both_file_obj.write(str(line))
        both_file_obj.close()
        for obj in write:
            print >>obj, '...%s total proteins written to %s' % \
                  (2*pcount, os.path.split(both_name)[1])
    #
    if MAKE_SEPARATE_FORWARD:
        for obj in write:
            print >>obj, '...%s proteins written to %s' % \
                  (pcount, os.path.split(for_name)[1])
    if MAKE_SEPARATE_REVERSED:
        for obj in write:
            print >>obj, '...%s proteins reversed and written to %s' % \
                  (pcount, os.path.split(rev_name)[1])
    #
    # close files and delete unwanted files
    #
    for_file_obj.close()
    rev_file_obj.close()
    fasta_lib.time_stamp_logfile('>>> ending: reverse_fasta.py', log_obj)
    log_obj.close()
    if not MAKE_SEPARATE_FORWARD:
        os.remove(for_name)
    if not MAKE_SEPARATE_REVERSED:
        os.remove(rev_name)
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
        database = r'E:\KUR1502\databases'
        fasta_file = fasta_lib.get_file(database, \
                                        [('FASTA files', '*.fasta'),\
                                         ('Zipped files', '*.gz'),\
                                         ('All files', '*.*')],\
                                        'Select a FASTA database')
        if fasta_file == '': sys.exit() # cancel button response
    #
    # call reverse function
    #
    os.chdir('.')
    print os.getcwd()
    fasta_reverse(fasta_file)
    try:    # wait for user to end program if not running from IDLE
        # what __file__ is under XP IDLE: 'C:\\Python26\\Lib\\idlelib\\idle.pyw'
        if not __file__.endswith('idle.pyw'):
            raw_input('\n...hit any key to end program...')
    except NameError:
        pass
#
# end
#
