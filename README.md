# fasta_utilities
## A collection of Python scripts for managing protein FASTA files

There can be many steps in getting a current FASTA database and preparing it for use by a search engine.  The database has to be downloaded to an appropriate location on your computer.  The database may need to be renamed to include version numbers.  The database will need to be uncompressed to a text file.  Many database releases include proteins from a large number of species and specific species subset databases (human, mouse, etc.) may need to be extracted.  Common contaminant sequences (trypsin, keratins, etc.) often need to be added.  Decoy sequences need to be created and appended to the databases when using that method of error estimation.  The version numbers, download dates, and number of protein sequences should be recorded and saved for inclusion in subsequent publications.  It is easy to make mistakes with so many steps or delay updating databases because it is too much work.  Hopefully, this set of utilities addresses many of these issues.

### Scripts and descriptions

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


## Scripts require Python 3

### Installing

There is no installation per se. The scripts use Python 3 and the only requirement is to have a reasonably recent Python 3 distribution installed on your system. You can use a [distribution from Python.org](https://www.python.org/downloads/), or you can install a scientific Python distribution like [Anaconda](https://www.anaconda.com/distribution/). These scripts do not use any functionality beyond what is in the standard distribution of Python. There are many other useful Python extensions in a scientific Python distribution like Anaconda and that is what I would recommend. If you are new to Python, you can start with a Python from Python.org and switch to Anaconda later. Installation and use of the standard distributions from Python.org may be a little simpler.  

Download this repo as a .zip file or clone it to access the script files. The scripts have been tested on Windows 7, Windows 10, and macOS 10.14. The scripts use some functions and class definitions in "fasta_lib.py", and some of the GUI scripts need "reverse_fasta.py" and a contaminants FASTA database. The scripts can be located anywhere on your computer. The fasta_lib script, other fasta utility scripts, and associated contaminant FASTA files should all be located inside of the same folder.

There are many ways to run Python scripts that depend a little on computer platforms and particular Python distributions. Python has a built-in development environment called IDLE that I like to use for most coding. The [IDLE documentation](https://docs.python.org/3/library/idle.html), a [shorter YouTube video](https://www.youtube.com/watch?v=Y47IzX4k6Go), or [longer video](https://www.youtube.com/watch?v=kXbpB5_ywDw) are avaialble for help. There are other IDEs that are popular, too. Anaconda distributes Sypder, for example. The Python distributions from Python.org can make it easier to find and use IDLE. One way to launch IDLE with Anaocnda is to open an Anaconda prompt window (an OS shell window with some environment variables defined) and type "idle3" on Mac or "idle" on Windows. Screen shots and how to run program instructions **below** will be using IDLE. Other ways to run Python and Python scripts are also possible. Experiment a little and find the method that you prefer.

IDLE has a console window that runs the main Python shell. This window can be resized, saves lines off of the top, and supports the clipboard (you can save output). IDLE includes a syntax highlighting editor that supports multiple windows. There is also GUI debugging functionality. Using IDLE to run the scripts goes something like this:

- launch IDLE
- open the file of the script you want to run
- run the script (run menu or F5)
- browse to FASTA files using the dialog boxes  

### Documentation

This README file is the main documentation for the scripts. The utilities were first written in 2010 using Python 2. The original documentation and a poster presented at the 2010 ASMS meeting are located in a "2010_documentation" folder. Much of the informations in the older documentation is still useful.

## NOTE

Several of the scripts access FTP sites at UniProt, NCBI, or Ensembl using Python standard library functions. Those are not very robust to erratic connections. They do not have built-in error handling and retry options. If you get error messages in the console output, a poor connection is the likely reason. The scripts have actually been tested and they tend to work most every time.

This is not professional software development (sorry). This is academic software with its typical limitations. I do not have a rigorous software testing suite. I do not have an update and release schedule. I fix bugs as I find them when using the software and push the changes. If you have some odd behavior with a script, try downloading a fresh copy. Please email me (wilmarth *AT* ohsu.edu) if you have any problems.

# User Guide

This guide describes a collection of Python scripts that automate downloading, organizing, naming, and managing FASTA protein databases.

## Background

The protein database choice is as important as your search engine choice and parameter settings. There are many major sources for protein databases (UniProt, NCBI, and Ensembl), various types of sequence collections (reviewed/processed sequences, reference proteomes, isoforms and variants, etc.), and more specialized collections (flybase, wormbase, yeast, tritryp, etc.). There can be pros and cons to specific protein databases that depend on the experimental goals. Identifying the proteins in your sample is usually just one step in the interpretation of the results. What you can do with the lists of identified protein can have a strong dependence on which protein database you used. You will need more than a passing knowledge of protein databases to do good proteomics research.

> Spectral libraries are often annotated with search engine identifications that will depend on the protein database used.  

The developers at Matrix Science have some nice descriptions of protein databases from several sources at [this link](http://www.matrixscience.com/help/seq_db_setup.html). There are some general descriptions of databases including details about FASTA header line formats.  This is a good time to talk about the general format for [protein FASTA files](https://en.wikipedia.org/wiki/FASTA_format). Here are a couple of FASTA entries from UniProt:

```
>tr|A0A087WNZ6|A0A087WNZ6_MOUSE RIKEN cDNA 1700025B11 gene OS=Mus musculus OX=10090 GN=1700025B11Rik PE=4 SV=1
MSTKNEEQNEDQSESVVIPHIQDHHCLAILAFCFFFPLGYLAYRSSCKTRTYIEQKEYEKAKATSRCTFA
FVLSSIAGGSIIFFCLFSRLFFM
>tr|A0A087WP11|A0A087WP11_MOUSE Histone H2A OS=Mus musculus OX=10090 GN=H2al1b PE=3 SV=1
MAKKMQRRRRQKRTRSQRGELPLSLVDRFLREEFHSSRLSSSALSFLTSVLEYLTSNILELAGEVAQTTG
RKRIAPEDVHLVVQNNEQLRQLFKPGGTSVNEDDN  
```

## FASTA entry format

FASTA files have header lines that start with a ">" character and they can be very long. After the header line is the protein sequence that is typically one or more lines long. These lines are usually of a fixed length (60, 70, or 80 characters are common), consist of valid amino acid characters, and are usually in upper case. There is no formal end-of-file designation. End of line characters can vary by operating system and can cause problems when FASTA files are moved between different platforms. The generic format of the header line is that the first character must be the ">" followed by a contiguous string of characters that function as a unique key for that entry in the file. This key string is often called the "accession". The rest of the line can be literally anything or even nothing. If there is additional information besides the accession, there needs to be a space (or some other form of whitespace) between the accession and the remainder of the line. In summary, a FASTA file consists of one or more FASTA entries. Each FASTA entry consists of a single header line (that starts with a ">" character) and an amino acid sequence (for proteins) that can span one or more lines. The header line consists of a contiguous string of characters (the accession that functions as a unique key) and an optional "description" string (with no defined format).

The general reading logic for a FASTA file goes something like this:
- open file and start reading lines
- if a line starts with ">" then a new (or the first) entry is starting
- strip the ">" character from the header lines
- find the first whitespace character
- characters up to the first whitespace character are the accession
- any characters after the first whitespace up to the end-of-line character are the description string
- read and save any lines until the next line that starts with ">" or end-of-file
- make the sequence string from the collection of lines
- save the sequence, accessions, and description in an appropriate data structure

After reading a FASTA file, there will be a collection of sequence entries. Each entry will have a protein sequence string, an accession key, and an (optional) description string. Each accession should be unique to serve as a potential key in a hashed data structure. The scope of uniqueness has to be at least the original FASTA file, but there may be reasons to extend the scope (e.g. a search engine that accepts multiple FASTA files). There are typically additional structure and/or restrictions for each of these parts of a FASTA entry.

### Protein Sequence

The amino acid sequence can consist of the usual 20 amino acid characters, and a few other characters. There are two less common amino acid characters (O and U), and some symbols for ambiguous situations. Some older sequencing methods convert amines to acids so that N and D are indistinguishable (as are Q and E). The symbols B and Z denote the two indistinguishable amide/acid cases, respectively. J denotes indistinguishable I/L. X denotes one (or more) unknown amino acids. The are two valid special characters: "\*" denotes a stop codon (these are not always at the end of sequences - and most sequences do not end in a \*), and "\-" denotes a gap of unknown length. One of the great mysteries in life is what does each of the popular search engines do when they encounter these valid characters that are not one of the 20 canonical residue symbols. B and Z can be handled by splitting the difference in the residue masses that differ by 0.984 Da. This worked okay when we did not have high mass accuracy instruments. I and L are the same mass so J is a pretty easy edge case. X is a good head scratcher because they are many possibilities in theory and (it seems) each search engine tries to pick a different one to implement. SEQUEST from Thermo uses Leucine in place of X. Mascot tries all 20 amino acids. I do not know what other search engines do with X. The situations is even more unclear for "\*" and "\-".

Different sources for protein databases can have different likelihoods of having these less common characters present in the amino acid sequences. One strategy is to test protein databases for these characters and see if there is anything to worry about. Filtering sequences before searching might be needed if the search engine has issues with certain characters. We will leave this topic unresolved and move on.

### Accessions

The accession is often a composite string composed of multiple elements separated by a joining character ("|" and "\_" are common). Some of these elements can be relatively unique and some not so much. In the above example sequences, [UniProt accession format](https://www.uniprot.org/help/fasta-headers) is a database designation, a stable accession, and a more human friendly identifier (e.g. 'tr|A0A087WNZ6|A0A087WNZ6_MOUSE'). The "tr" denotes a TrEMBL entry (an unreviewed entry; "sp" denotes a reviewed Swiss-Prot entry), "A0A087WNZ6" is the stable accession string, and "A0A087WNZ6_MOUSE" is the concatenated gene symbol and species OS code. NCBI and other database sources use different accession formats that are described in their manuals and help pages. Decoy sequences often have modified accessions based on the original accessions with prefix or suffix elements. Like most composite database keys, the constituent elements carry no guarantee of uniqueness. Only the combination of the constituent elements is guaranteed to be unique. Historically, some sort of regular expression parsing of accessions seems to be the norm. This can have consequences because the simplified parts of the accession that are retained may not be unique. If parsing is done before search engine processing, the parsed accessions need to be tested for uniqueness. It is also common in database searching to add common laboratory contaminant sequences and decoy sequences for error estimates. Any regular expression parsing rule has to be flexible enough for all types of accessions that will be encountered (contaminants can be the most variable), or the extra FASTA entries (contaminants and decoys) have to be formatted to have accessions compatible with the parsing expression. It is far safer to not parse accession strings until after search engine results have been compiled into protein results lists.

### The rest of the header line

If the sequences and accessions are not complicated enough to keep you awake at night, do not worry. The description strings more than make up for it. They are have no defined format and are only limed by the imagination of database creators, as illustrated by the [PEFF](http://www.psidev.info/peff) proposal. UniProt includes several bits of information beyond the actual protein/gene description phrase. There is commonly the species (OS=), the taxonomy number (OX=), the gene symbol (GN=), the protein evidence code (PE=), and version number (SV=). NCBI is know for its non-redundant protein database. Those databases have one sequence for all species where that protein is the same (e.g. 138 higher eukaryotic ubiquitin sequences are all identical). The NCBI entry contains all 138 FASTA header lines in one composite header line separated by Control-A characters. I have seen FASTA header line lengths of several thousand characters. Ensembl protein databases have complicated descriptions that carry a full compliment of cross-reference information on corresponding genome and transcriptome coordinates, in addition to gene symbols, protein descriptions, and sequence sources and types. The accessions also have a version number suffix.  An example sequence for human ENSP00000374886 is shown below from Ensembl release 93. Other protein database sources may have other types of information encoded into their FASTA header lines. Some detailed research is required for each source to see what information is present and if it is in a usable format for a protein researcher. There can be benefits to reformatting FASTA header descriptions into more concise and relevant strings.

```
>ENSP00000374886.2 pep chromosome:GRCh38:7:142391891:142392412:1 gene:ENSG00000211716.2 transcript:ENST00000390363.2 gene_biotype:TR_V_gene transcript_biotype:TR_V_gene gene_symbol:TRBV9 description:T cell receptor beta variable 9 [Source:HGNC Symbol;Acc:HGNC:12246]
MGFRLLCCVAFCLLGAGPVDSGVTQTPKHLITATGQRVTLRCSPRSGDLSVYWYQQSLDQ
GLQFLIHYYNGEERAKGNILERFSAQQFPDLHSELNLSSLELGDSALYFCASSV
```     

## Protein database sources
The best sources for protein sequence collections depends a little on whether you are studying model systems, eukaryotes, prokaryotes, plants, parasites, etc. Some sources are genomic oriented, some are protein oriented. Some sources are dedicated to specific systems and serve as centralized resources for researchers. There are many sources of protein databases and all of them seem to be in use to some degree. The major sequence collections are described in this [Wikipedia entry](https://en.wikipedia.org/wiki/Sequence_database). [UniProt](https://www.uniprot.org/) databases are probably the most widely used protein databases in proteomics with databases from [NCBI](https://www.ncbi.nlm.nih.gov/protein) the next most common. [Ensembl](https://uswest.ensembl.org/index.html) databases can be useful if protein information is being compared with genome and transcriptome information. Ensembl genomes are often references for read alignments in next generation sequencing.

There are several sources for commonly studied model systems such as [yeast](https://www.yeastgenome.org), [fruit flies](http://flybase.org), [worms](https://www.wormbase.org/#012-34-5), [disease vectors](https://www.vectorbase.org), [trypanosomes](http://tritrypdb.org/tritrypdb/), [zebrafish](https://zfin.org), and [Arabidopsis](https://www.arabidopsis.org/index.jsp). Many times these dedicated resources are the sources of the sequence collections used in the more centralized database sources. Dedicated sources provide FASTA protein database download options in many cases. It is wise to do some investigation of these databases for accession and description formats, and to see if unusual amino acid sequence characters are an issue. The centralized sources of protein databases (UniProt, NCBI, and Ensembl) have processing pipelines so that FASTA header lines have defined and consistent formats, and have had some consistent set of rules applied to the amino acid sequences (e.g. how are stop codons handled?). As the Grail Knight in [Raiders of the Lost Arc](https://www.imdb.com/title/tt0082971/) might say, "choose your protein database wisely."

Protein databases from different sources (and different databases options from the same source) can be highly variable in both total protein sequence count and the degree of peptide redundancy (we usually work with tryptic peptides). For the curious, this [Master's thesis](https://digitalcommons.ohsu.edu/etd/3855/) by Ravi Madhira explores this topic in some detail.

## UniProt sequence collections
[UniProt](https://www.uniprot.org/help/about) has a lot of protein sequences that are grouped and organized in a variety of ways. The sequences are composed of [two sections](https://www.uniprot.org/help/uniprotkb_sections): manually reviewed Swiss-Prot sequences and computer annotated TrEMBL (unreviewed) sequences. These are denoted with "sp" and "tr" prefixes, respectively, in FASTA file accessions. UniProt has a monthly release schedule for its databases.

UniProt pulls sequencing data into its database from other sources daily. The sequences are run through a series of processing algorithms to produce the TrEMBL entries. Sequences in UniProt are organized by species and are non-redundant by species. New sequences are associated with their respective species and added to the sequence collections if they are unique. If there are protein variants that have differing sequences, these will appear as different entries in TrEMBL. The number of proteins being sequenced greatly exceeds that capacity for human reviewers to add them to Swiss-Prot, so most proteins in the UniProt collection are TrEMBL entries. Swiss-Prot entries are manually curated and they are fundamentally different from TrEMBL entries. Some subset of the annotation information for each protein can be obtained programmatically and that information tends to be present in TrEMBL and in Swiss-Prot entries. There is much more information added to the Swiss-Prot database records during the curation process.

It is important to understand the relationship between Swiss-Prot and TrEMBL databases at UniProt. They are not redundant; they are orthogonal. Information for a specific protein exists in either Swiss-Prot or in TrEMBL, but not in both. Although TrEMBL is much larger than Swiss-Prot, Swiss-Prot is the tail that wags the dog. It is more correct to think of TrEMBL as a queue that holds proteins until they can be curated and added to Swiss-Prot. The rate that entries migrate from TrEMBL to Swiss-Prot rivals the rate that wealth trickles down after tax breaks to the 1%, but the historical relationship between Swiss-Prot and TrEMBL is what we have at present. When information is added to Swiss-Prot, information disappears from TrEMBL.

Protein variants are formally annotated as deltas relative to a chosen canonical sequence. The canonical sequence is not chosen based on biological relevance; it is chosen for annotation efficiency. The longest sequence is often used because describing shorter forms is easier. Because variants are described relative to a base sequence, several TrEMBL entries can get rolled into one Swiss-Prot record. For a given species, the protein sequences actually fall into three distinct categories: there are Swiss-Prot canonical sequences (one protein sequence for each Swiss-Prot record that is associated with one accession), there are isoforms/variants of those Swiss-Prot proteins, and there are almost always additional TrEMBL entries awaiting curation.

UniProt gives you options to download sequences in many different ways. You can access the main web interface and search for proteins. Any matching proteins can be further filtered by options like "reviewed" or "unreviewed" sequences. The matching, filtered proteins can be downloaded in FASTA format. There may be options to download canonical sequences or canonical and isoform sequences. There are also predefined sequence collections of various sorts (proteomes, reference proteomes, etc.) available via the web interface. There is also a UniProt FTP site where FASTA files can be downloaded. There can be protein databases available via FTP that do not exist via the web interface. Full (all species) sequence collections is one example. All Swiss-Prot or all UniProt sequences are available. Swiss-Prot is around 600K sequences. The full sequence collection is very much larger (a 32 GB compressed file with over 133 million sequences as of 11/10/2018 from close to one million species). There are scripts here that can download and analyze these large sequence collections. There are other scripts that can extract FASTA protein databases from these collections by species or groups of species. There is also a GUI script that downloads reference proteomes from the FTP site.

## National Center for Biotechnology Information (NCBI) sequences
[NCBI](https://www.ncbi.nlm.nih.gov/) is a major repository for genomic and protein sequences. Being a collector of every sequence under the sun poses some logistic challenges. Without some cleanup of the sequences, the sheer volume of information is somewhat overwhelming. The solution to this computer generated problem is, of course, computer automated processing via the [RefSeq project](https://www.ncbi.nlm.nih.gov/refseq/). Pretty much any database from NCBI that you would want to use with a search engine on proteomics data should be a RefSeq database. NCBI also has protein databases organized in various sets of related sequences (e.g. human reference genome, prokaryotic RefSeq genomes, etc.) that can be retrieved via different web site links and an extensive FTP site. I think the sequence collections at NCBI are continuously updated, i.e. no formal versioning for the whole kit and caboodle. Individual sequences have version number suffixes on accessions.

NCBI is famous for the BLAST algorithm and that is powered by the infamous NCBI nr protein database. The nr stands for non-redundant. This database, which can be downloaded from the FTP site, is basically one of every protein sequence currently known to man (and other genders). The FASTA header lines are epic, compound lines (multiple headers separated by Control-A characters) to carry complete annotations for the proteins. Extracting sequences by species involves testing each part of the FASTA header lines to see if the species of interest has been associated with the respective sequence. Interestingly, when downloading protein databases from web browsers, the databases are **not** non-redundant. NCBI accessions also have a [complicated format](https://www.ncbi.nlm.nih.gov/Sequin/acc.html) that covers genomic and protein accessions. There is a script that will fetch the nr database from the FTP site and analyze the sequences by species, and a couple of scripts for extracting sequences. The nr database is large (42 GB compressed on 11/11/2018) as is the associated species mapping information (593 million species identified accessions). The compressed FASTA file has over 177 million sequences (284 million for non-redundant by species which would be similar to how UniProt counts things) from 1.15 million different species. It is no longer possible to keep all of the necessary data structures in memory, so the analysis of the nr FASTA file to get sequence counts by species is now much slower. The continued growth of the NCBI nr database makes this method of getting sequences obsolete (lots of memory and SSD drives can help). At this point, the logical strategy would be to explore other download mechanisms from NCBI. There are no GUI applications for getting single species protein databases from the FTP site (like the UniProt and Ensembl tools) yet. This would be a great summer student project.

NCBI made a [major change](https://ncbiinsights.ncbi.nlm.nih.gov/2016/07/15/ncbi-is-phasing-out-sequence-gis-heres-what-you-need-to-know/) to its sequence identifiers in 2016. For years (decades?) the primary key was a "gi" number (GenInfo Identifier). The gi numbers have been replaced by compound accession "dot" version [format](https://ncbiinsights.ncbi.nlm.nih.gov/2016/12/06/converting-gi-numbers-to-accession-version/). If you use any of these major database sources, formats and options do change from time to time, and it pays to visit website documentation frequently to see if there are changes.

## Scripts for downloading and extracting databases
These scripts get protein FASTA files from respective sites:

- nr_get_analyze.py
- sprot_get_analyze.py
- uniprot_get_analyze.py

They all import the fasta_lib.py module (a collection of common functions and classes). They fetch large multi-species databases from FTP sites. The NCBI nr data has grown so large that the nr_get_analyze script is painful to run. The size of nr is still growing exponentially. More direct ways to get to relevant subsets of nr are available, although it would take some effort to wrapper the choices to make finding, organizing, and naming downloads convenient. The Swiss-Prot database from UniProt has not grown so quickly and extracting species from this database is still manageable. TrEMBL is large, but UniProt has done more to control the size growth. Working with Siwss-Prot and TrEMBL is still possible but it will stress test your computer and internet connection.

After the above scripts have downloaded their respective databases, the number of sequences for each species are tallied and written to tab-delimited files that can be opened with a spreadsheet program to provide the necessary information to decide what sequences to extract for database searching. A support script "taxon_group_analyzer.py" provides summaries of the sequences associated with taxonomy "nodes" (e.g. rodent).

There are companion extraction scripts for the multi-species database downloads listed above that can extract by taxonomy number or by text strings:

- extract_by_string.py
- nr_extract_taxon.py
- uniprot_extract_from_one.py
- uniprot_extract_from_both.py

There are two UniProt scripts because it makes sense to get just Swiss-Prot sequences (for some species) **or** sequences from **both** Swiss-Prot and TrEMBL. You **never, ever** want to use just TrEMBL sequences. Sequences of any proteins present in Swiss-Prot for a respective species will have been removed from TrEMBL during curation. **Note:** any scripts that start with lowercase were part of the original 2010 utilities and the Word files in the 2010_documentation folder will have more detailed documentation.

In 2017, a summer student (Delan Huang) and I created a couple of GUI scripts to help get [reference proteomes](https://www.uniprot.org/help/reference_proteome) from UniProt and to get [Ensembl vertebrate](https://uswest.ensembl.org/index.html) databases.

- UniProt_reference_proteome_manager.py
  - fasta_lib.py
  - reverse_fasta.py
- Ensembl_proteome_manager.py
  - fasta_lib.py
  - Ensembl_fixer.py

## Ensembl Proteome Manager

- Ensembl_proteome_manager.py
- fasta_lib.py
- Ensembl_current_release.pickle
- Ensembl_fixer.py
- default_Ensembl_species.txt

This script uses a GUI window (in addition to some console output) to show you the list of vertebrate species in the current Ensembl release (149 proteomes as of 11/12/2018). There are options to filter the list of proteomes to find those of interest. And options to add contaminants and/or decoys to the downloaded databases. Different contaminant databases can be used. The list of downloaded proteomes can be saved so that those species can be updated more easily. File and folder naming is done automatically to append release information and keep the downloaded databases organized. The FASTA header lines in Ensembl databases are not very friendly for typical researches (my opinion) and they are reformatted and shortened to be more useful.

![Ensembl Main GUI window](/images/Ensembl_1_main_edited.jpeg)

**Ensembl main window.** The GUI has a frame at the top to facilitate searching for proteomes and for specifying how to process the downloaded FASTA files. The lower frame has a left side with the available proteomes and a right side with the desired proteomes to download. The center set of buttons manage the left and right lists and the downloading. There is a status bar at the bottom.

---

![Ensembl filtering controls](/images/Ensembl_2_top_edited.jpeg)

**Filtering the proteome list and processing options.** The left list of proteomes can be filtered based on species names or taxonomy numbers. The searching is not case sensitive and does a simple "in" test. Substrings will return results and the general idea is to make the left list a little shorter so the species of interest can be found more easily. The downloaded databases will be compressed. During decompression, common contaminants can be added from a specified contaminants FASTA file. Sequence reversed decoys can also be added. The two check box options are independent and both can be checked.

---

![Ensembl filter for mouse](/images/Ensembl_3_mouse_edited.jpeg)

**Example of how to get mouse proteomes.** The taxonomy number for mouse is 10090. If we enter that in the TAxonomy ID field, and click the Show Filtered List button, we will get 13 mouse proteomes. Ensembl has specific proteomes for 13 common mouse strains. The top one in the left list is the typical mouse genome of the most commonly used strain (C57BL/6J, I think).

---

![Ensembl selecting downloads](/images/Ensembl_4_add_edited.jpeg)

**Adding mouse to the download list with human.** If we select the first mouse line on the left, then click the Add Proteome(s) button, that proteome is added to the right window. We can click the Download button to download and process these databases.

---

![Ensembl download dialog](/images/Ensembl_5_download_edited.jpeg)

**A dialog box lets you select the location for Ensembl databases on your computer.** The script will take care of creating release version named subfolders. You want to pick a "container" folder for your collection of Ensembl databases. Examination of the subfolders and their contents will give you an idea of the general organization and naming scheme. Some information in the filenames is redundant with information in the folder names on purpose. When adding FASTA files to data repositories or as Supplemental files, all of the release information should be contained in the filename because the file is usually taken out of it folder path context.

---

![Ensembl file organization](/images/Ensembl_6_files.png)

**Subfolder organization.** The compressed downloaded files from the Ensembl FTP site are located in the folders with species information. The decompressed databases have the ".fasta" file extensions. We did not select any processing options, so we just have the target databases without any common contaminants. There is also a log file with the information that was shown in the console window when the script ran. This window also has one of the mouse strains (left over from an earlier testing).

## UniProt Reference Proteome Manager

- UniProt_reference_proteome_manager.py
- fasta_lib.py
- UniProt_current_release.pickle
- default_UniProt_species.txt
- reverse_fasta.py

UniProt has several ways to find and download databases. The main web site options are the easiest to find and use. They have limitations, however. There is a [UniProt FTP](https://www.uniprot.org/downloads) site that is often overlooked. There is a reduced list of higher quality reference proteomes, for example. There are (as of 11/15/2018) 439 archaea, 8895 bacteria, 1184 eukaryota, and 6178 virus reference proteomes available via FTP. The sequence collections for each species are split into a canonical set (sort of a one gene one protein idea) and (optionally) any additional isoforms of canonical proteins.

Protein databases available through the main web site are split into reviewed sequences (Swiss-Prot entries) and unreviewed entries (TrEMBL entries). Swiss-Prot (and only Swiss-Prot) entries can have optional annotated isoforms. The canonical sequence collections contain both Swiss-Prot and TrEMBL entries to make up "complete" proteomes. The higher eukaryotic canonical proteomes all have around 21000 sequences, for example. These canonical databases are, therefore, reasonably complete with minimal peptide redundancy. These databases are particularly good choices for shotgun quantitative proteomics data.

There are README files that provide the mappings from the species names and taxonomy numbers to the UniProt proteome numbers. The actual FTP file listings only have the proteome numbers. This can make finding the right databases to download a little challenging. This script gets information from the FTP site and presents it in a more human friendly format. There is automatic folder creation and file naming logic to help keep databases organized and make sure that the UniProt release information is captured. There are also some convenience options to add common contaminants and decoy sequences. _**Note:** The folder naming options are not yet implemented._      

![UniProt main window](/images/UniProt_1_main_edited.jpeg)

**UniProt reference proteome manager main window.** There is a top pane for controlling what proteomes are presented in the left list of the lower pane. Different kingdoms can be selected, species names can be restricted to specific species names, and taxonomy numbers can also be restricted to those of interest. After downloading databases, they can be processed to add contaminants or decoys (and contaminants). The user can select different contaminant databases if desired.

The bottom pane has available proteomes listed on the left, and the desired databases to download on the right. The right list can be saved as a default list that loads when the GUI launches. There are controls to move proteomes from the left to the right, to drop proteomes from the right list, and download the databases. Databases can be downloaded as canonical sequences only, or canonical sequences and isoform sequences. If isoforms are downloaded, they will be automatically added to the canonical sequences to make a single combined protein FASTA database.

---

![UniProt bovine species search](/images/UniProt_2A_bovine_edited.jpeg)

**Filtering for bovine (cow) sequences.** We can restrict the kingdom to Eukaryota by unchecking the other kingdom boxes. We can require that the species name contain "bovine" (a case insensitive "in" test), and click the Show Filtered List button.

---

![UniProt left list filtered](/images/UniProt_2B_bovine_edited.jpeg)

**Left list will update.** We now have just two possible proteomes on the left. We can select the bovine proteome (taxonomy 9913) and click the Add Proteome(s) button.

---

![UniProt moved to right](/images/UniProt_2C_bovine_edited.jpeg)

**Bovine proteome added to right list.** The bovine proteome has been added to our download list. The right list has some species loaded from our default list that we do not need.

---

![UniProt select to drop](/images/UniProt_3A_drop_edited.jpeg)

**Drop any unneeded proteomes.** We can select the mouse, rat, yeast, and E. coli rows and then click the Drop Proteome(s) button to remove them.

---

![UniProt right updated and download](/images/UniProt_3B_drop_edited.jpeg)

**Ready to download databases.** We are ready to download some protein databases to test if there are human proteins that make us behave more like a herd of cattle, a flock of sheep, or if we really are just a bunch of dirty little pigs. We will download just the canonical sequences and we will add decoys and contaminants so the databases are ready to use with a search engine program, such as [Comet](http://comet-ms.sourceforge.net/).

---

![UniProt download dialog](/images/UniProt_4A_download_edited.jpeg)

**Specify the download location.** We want to select a folder where we will keep all of our UniProt protein database. Subfolder creation, naming, and file naming will be taken care of by the script.

---

![UniProt console](/images/UniProt_4B_download_edited.jpeg)

**Console window also has information.** Download progress is logged to the console window with details on filenames, locations, and sequence counts.

---

![UniProt after download and quit](/images/UniProt_4C_download_edited.jpeg)

**Quit after downloads finish.** The status bar and a dialog alert box will let you know when downloads have finished. If you do not have any additional databases to download, it is time to quit. Click the quit button or close the GUI window. You may also want to quit your Python 3 shell.

---

![UniProt update defaults](/images/UniProt_5_defaults_edited.jpeg)

**Right list can be saved.** The right list might be a collection of species that you want to download on a regular basis. If the current right list differs from the previously saved list, you will be asked if you want to save the changes.

---

![UniProt files](/images/UniProt_6_files_edited.jpeg)

**Example of what downloaded files/folder look like.** The compressed download files are saved in nicely named folders. The downloaded files are decompressed, descriptively named, and any selected processing performed. A FASTA file of the downloaded database is always created in addition to any desired processed versions (with decoys or contaminants). A log file is also present.

## Scripts for Working with Downloaded FASTA Files

### add_extras_and_reverse.py
Adds extra sequences to protein databases. The extra sequences need to be in a separate FASTA file (often only a few sequences). The accessions of the extra proteins are modified to avoid any accession conflicts. Contaminants and decoys can be added to the resulting FASTA file.

### check_for_duplicates.py

Checks a FASTA file for duplicated protein sequences. Produces a report of any duplicates that were found.

### count_deluxe_fasta.py

Counts the sequences in one or more FASTA files with valid amino acid character testing. Produces a report with sequence lengths and calculated molecular weights.

### count_fasta.py

Counts the sequences in one or more FASTA files.

### Ensembl_fixer.py

Reformats FASTA header lines in Ensembl protein databases into a more human-readable, concise line.

### FASTA_digester.py

Performs a theoretical digest of a FASTA protein database and produces a report of peptide redundancy (and some other statistics). The default is a tryptic digest. Script modification is necessary to support other proteases. The set of digestion options for [Comet](http://comet-ms.sourceforge.net/) are available.

### remove_duplicates.py

Creates a non-redundant protein database along the same lines as the nr release from NCBI. Duplicated sequences will appear once with a compound FASTA header line separated by Control-A characters.

### reverse_fasta.py

Adds reversed decoy sequences (and contaminants) to FASTA files. Concatenated (recommended) or separate decoy database can be produced.

### TriTryp_fixer.py

Does some protein sequence character checking and FASTA header line reformatting of protein data bases from [TriTryp](http://tritrypdb.org/tritrypdb/).

---
#### Details
This documentation written by Phil Wilmarth, OHSU, November 2018.
