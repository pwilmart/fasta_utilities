# fasta_utilities
## A collection of Python 3 scripts for managing protein FASTA files.
### Winter 2018

- Ensembl_fixer.py - does header line reformatting for v83 and newer Ensembl fasta databases
- Ensembl_proteome_manager.py - GUI for downloading Ensembl fasta databases
- FASTA_digester.py - theoretical digestion statistics of protein databases
- TriTryp_fixer_v2.py - reformats fasta header lines and does some sequence analysis
- UniProt_reference_proteome_manager.py - GUI for downloading UniProt reference proteomes from the FTP site
- add_extras_and_reverse.py - adds special sequences to protein databases
- check_for_duplicates.py - looks for duplicate protein sequences
- count_deluxe_fasta.py - counts entries in multiple FASTA files
- count_fasta.py - counts entries in FASTA file
- extract_by_string.py - creates subset databases by header line string patterns
- fasta_lib.py - main library module
- nr_extract_taxon.py - extracts subset databases from NCBI nr by taxonomy numbers
- nr_get_analyze.py - downloads and analyzes NCBI nr releases
- remove_duplicates.py - removes duplicate FASTA entries
- reverse_fasta.py - does simple sequence reversal for decoy generation
- sprot_get_analyze.py - downloads and analyzes UniProt Swiss-Prot releases
- taxon_group_analyzer.py - analyzes databases by taxonomy node numbers
- uniprot_extract_from_both.py - extracts species by taxonomy number from Swiss-Prot + TrEMBL databases
- uniprot_extract_from_one.py - extracts species by taxonomy number from Swiss-Prot databases
- uniprot_get_analyze.py - downloads and analyzes UniProt Swiss-Prot + TrEMBL releases


## Utility scripts have been re-written for Python 3.
The documentation (*.doc) files are out-of-date and have references to Python 2.7. An updated set of documentation is on the to-do list. Many utilities are the same (and function similarly) to what is described in the documentation. There are some new utilities and some new functionality to support newer options available at sources like UniProt and to handle changes in FASTA description/accession formats (NCBI and Ensembl).

# NOTE
Several of the scripts access FTP sites at UniProt, NCBI, or Ensembl using Python standard library functions. Those are not very robust to erratic connections. They do not have built-in error handling and retry options. If you get error messages in the console output, a poor connection is the likely reason. The scripts have actually been tested and they tend to work most every time. Please email me (wilmarth _AT_ ohse.edu) if you have any problems.
