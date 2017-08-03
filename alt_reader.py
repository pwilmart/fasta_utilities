import os
import gzip
import time
from itertools import groupby

# bare bones FASTA entry container
class Protein:
    def __init__(self):
        self.acc = ''
        self.desc = ''
        self.seq = ''

def fasta_iter(fasta_name):
    """Yields Protein objects from fasta files.
    Adapted from "https://www.biostars.org/p/710/" post.
    See also: https://drj11.wordpress.com/2010/02/22/python-getting-fasta-with-itertools-groupby/
    """
    with open(fasta_name) as fasta_handle:
        # skip x[0] (boolean something from groupby) and keep the alternating header, sequences
        # more on groupby here: https://docs.python.org/3/library/itertools.html
        fasta_iter = (x[1] for x in groupby(fasta_handle, lambda line: line[0] == '>'))
        for header in fasta_iter:
            p = Protein()
            header = header.__next__()[1:].rstrip()
            p.acc = header.split()[0]
            p.desc = header[len(p.acc)+1:]
            # join all sequence lines
            p.seq = "".join(s.strip() for s in fasta_iter.__next__())
            yield p
            
def fasta_gz_iter(fasta_name):
    """Yields Protein objects from fasta files.
    Adapted from "https://www.biostars.org/p/710/" post.
    See also: https://drj11.wordpress.com/2010/02/22/python-getting-fasta-with-itertools-groupby/
    """
    with gzip.open(fasta_name, mode='rt') as fasta_handle:
        # skip x[0] (boolean something from groupby) and keep the alternating header, sequences
        # more on groupby here: https://docs.python.org/3/library/itertools.html
        fasta_iter = (x[1] for x in groupby(fasta_handle, lambda line: line[0] == '>'))
        for header in fasta_iter:
            p = Protein()
            header = header.__next__()[1:].rstrip()
            p.acc = header.split()[0]
            p.desc = header[len(p.acc)+1:]
            # join all sequence lines
            p.seq = "".join(s.strip() for s in fasta_iter.__next__())
            yield p

fasta_name = '/Users/pwilmart/Desktop/TimeMachineSkip/Databases_Spring2017/uniprot_2017.07/uniprot_trembl_2017.07.fasta.gz'

print('start', time.ctime())
proteins = 0
for protein in fasta_gz_iter(fasta_name):
    proteins += 1

print(proteins, 'proteins were read')
print('end', time.ctime())
    
"""Better way to have the function take a file handle instead of name. Then gz test is in the caller."""
