"""'taxon_group_analyzer.py' Written by Phil Wilmarth, OHSU.

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
import tarfile
import fasta_lib_Py3 as fasta_lib


def main(node_taxon):
    """Program to process taxonomy nodes file and find groups of species.
    """
    print('=======================================================================')
    print(' taxon_group_analyzer.py, v1.1.0, written by Phil Wilmarth, OHSU, 2017 ')
    print('=======================================================================')
    
    # get the name of the database analysis text file
    default = r'C:\Xcalibur\database'
    if not os.path.exists(default):
        default = os.getcwd()
    ext_list = [('Text files', '*.txt'), ('All files', '*.*')]
    analysis_file = fasta_lib.get_file(default, ext_list, 'Select a species analysis file')
    if analysis_file == '': sys.exit() # cancel button response
    
    analysis_folder, short_file = os.path.split(analysis_file)    
    print('...making taxonomy nodes dictionary...')

    # may need to check if this works in Python 3
    archive_name = os.path.join(analysis_folder, 'taxdump.tar.gz')
    archive = tarfile.open(archive_name)
    nodes = archive.extractfile('nodes.dmp')
    
    # read file and save taxon to parent taxon mappings
    taxon_to_parent = {}
    while True:
        line = nodes.readline()
        line = line.decode('utf-8')
        line = line.rstrip()
        if not line:
            break
        else:
            line = line.rstrip()
        item = line.split('\t|\t')
        taxon_to_parent[int(item[0])] = int(item[1])    
    nodes.close()
    
    # open the fasta_analysis.txt file and find group members
    print('...scanning %s file...' % (short_file,))
    fasta_analyze = open(analysis_file, 'r')
    out_name = analysis_file.replace('.txt', '_' + str(node_taxon) + '.txt')
    out_file = open(out_name, 'w')
    line = fasta_analyze.readline().rstrip()
    print('Analysis of node:', node_taxon, file=out_file)
    line = line.replace('A2:', 'A3:')
    print(line, file=out_file)
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
                break
        tree.append(1) # add last lineage item
        if node_taxon in tree:  # see if desired node is anywhere in the list
            member += 1
            print(line, file=out_file)  # write lines of node members
    #
    fasta_analyze.close()
    out_file.close()
    print('...taxonomy node %s had %s members...' % (node_taxon, member))
    return


# check for command line launch and see if a taxonomy number was passed
if __name__ == '__main__':
    
    # if arguments make sure it is an integer
    node_taxon = 0
    if len(sys.argv) > 1:
        try:
            node_taxon = int(sys.argv[1])
        except:
            pass
    
    # if no argument or improper argument, ask for taxon number
    if not node_taxon:
        node_taxon = input('... enter a taxon node ID number > ')
        node_taxon = int(node_taxon)
    main(node_taxon)

# end
