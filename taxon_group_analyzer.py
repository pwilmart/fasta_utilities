"""'taxon_group_analyzer.py' Written by Phil Wilmarth, OHSU.
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
def main(node_taxon):
#====================
    """Program to process taxonomy nodes file and find groups of species.
    """
    import os
    import fasta_lib
    import tarfile
    #
    print '======================================================================'
    print ' taxon_group_analyzer.py, v1.0, written by Phil Wilmarth, OHSU, 2009.'
    print '======================================================================'
    #
    # get the name of the database analysis text file
    #
    default = r'C:\Xcalibur\database'
    if not os.path.exists(default):
        default = os.getcwd()
    ext_list = [('Text files', '*.txt'), ('All files', '*.*')]
    analysis_file = fasta_lib.get_file(default, ext_list, \
                                       'Select a species analysis file')
    if analysis_file == '': sys.exit() # cancel button response
    analysis_folder, short_file = os.path.split(analysis_file)
    #
    print '...making taxonomy nodes dictionary...'
    archive_name = os.path.join(analysis_folder, 'taxdump.tar.gz')
    archive = tarfile.open(archive_name)
    nodes = archive.extractfile('nodes.dmp')
    #
    # read file and save taxon to parent taxon mappings
    #
    taxon_to_parent = {}
    while True:
        line = nodes.readline()
        if not line:
            break
        else:
            line = line.rstrip()
        item = line.split('\t|\t')
        taxon_to_parent[int(item[0])] = int(item[1])    
    nodes.close()
    #
    # open the fasta_analysis.txt file and find group members
    #
    print '...scanning %s file...' % (short_file,)
    fasta_analyze = open(analysis_file, 'r')
    out_name = analysis_file.replace('.txt', '_'+str(node_taxon)+'.txt')
    out_file = open(out_name, 'w')
    line = fasta_analyze.readline().rstrip()
    print >>out_file, 'Analysis of node:', node_taxon
    line = line.replace('A2:', 'A3:')
    print >>out_file, line
    member = 0
    while True:
        line = fasta_analyze.readline() # read analyze text file line
        if not line:
            break
        else:
            line = line.rstrip()    
        tree = [] # list of taxon number lineage
        parent = line.split('\t')[1]
        try:
            parent = int(parent)
        except:
            continue
        while parent != 1:  # all lineages end with taxon=1
            tree.append(parent)
            try:
                parent = taxon_to_parent[parent]
            except KeyError:
                #print '   WARNING:', line
                break
        tree.append(1) # add last lineage item
        if node_taxon in tree:  # see if desired node is anywhere in the list
            member += 1
            print >>out_file, line  # write lines of node members
    #
    fasta_analyze.close()
    out_file.close()
    print '...taxonomy node %s had %s members...' % (node_taxon, member)
    return
#
# check for command line launch and see if a taxonomy number was passed
#
import sys
if __name__ == '__main__':
    #
    # if arguments make sure it is an integer
    #
    node_taxon = 0
    if len(sys.argv) > 1:
        try:
            node_taxon = int(sys.argv[1])
        except:
            pass
    #
    # if no argument or improper argument, ask for taxon number
    #
    if not node_taxon:
        node_taxon = raw_input('... enter a taxon node ID number > ')
        node_taxon = int(node_taxon)
    main(node_taxon)
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
