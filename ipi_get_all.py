"""'ipi_get_all.py' Written by Phil Wilmarth, OHSU.
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
#=======================================
def main(DB, folder, versions, entries):
#=======================================
    """Fetches current versions of all IPI databases.

    Arguments: DB is a list of the desired IPI databases,
    "folder" is full path name to IPI DBs.  "versions" is a dictionary
    of version numbers, "entires" dictionary has sequence counts.
    """
    import os
    import fasta_lib
    #
    print '=============================================================='
    print ' ipi_get_all.py, v1.01, written by Phil Wilmarth, OHSU, 2009.'
    print '=============================================================='
    #
    # create a log file to mirror screen output
    #
    log_obj = open(os.path.join(folder, 'fasta_utilities.log'), 'a')
    write = [None, log_obj]
    fasta_lib.time_stamp_logfile('\n>>> starting: ipi_get_all.py', log_obj)
    #
    # make sure the files are present or download if not
    #
    for db in DB:
        fasta_lib.download_ipi(db, folder, versions)
    #
    # print out some stats and exit
    #
    for db in DB:
        new_db = 'ipi.%s.%s.fasta.gz' % (db.upper(), versions[db],)
        for obj in write:
            print >>obj, '...%s contained %s entries...' % (new_db, entries[db])
    #
    fasta_lib.time_stamp_logfile('>>> ending: ipi_get_all.py', log_obj)
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
    import fasta_lib
    #
    # get path to IPI databases and call main function to download, etc.
    # check if folder path is passed on command line
    #
    versions, entries = fasta_lib.get_ipi_versions()  # get current versions
    if len(sys.argv) > 1 and os.path.exists(sys.argv[1]):
        container = sys.argv[1]
    else:
        #
        # browse to a container folder for IPI downloads
        # a subfolder will be created with name "ipi_humanversion"
        #
        database = r'C:\Xcalibur\database'
        if not os.path.exists(database):
            database = os.getcwd()
        container = fasta_lib.get_folder(database, \
                                       'Select folder for IPI downloads')
        if container == '': sys.exit() # cancel button response
    folder = os.path.split(container)[1].lower() # ignore case
    if folder.startswith('ipi_'):   # subfolder selected, move one level up
        container = os.path.split(container)[0]
    #
    # create ipi download folder with the version from HUMAN db
    #
    folder = os.path.join(container, 'ipi_' + versions['HUMAN'])
    if not os.path.exists(folder):
        os.mkdir(folder)
    #
    # get list of databases and call MAIN
    #
    db_list = versions.keys()
    main(db_list, folder, versions, entries)
    try:    # wait for user to end program if not running from IDLE
        # what __file__ is under XP IDLE: 'C:\\Python26\\Lib\\idlelib\\idle.pyw'
        if not __file__.endswith('idle.pyw'):
            raw_input('\n...hit any key to end program...')
    except NameError:
        pass
#
# end
#


                
    


