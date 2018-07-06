"""'fasta_lib.py' Written by Phil Wilmarth, OHSU.

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
# 6/5/2009 PW - minor bug fixes for processing nr databases
#             - changed how reversed accessions are written for nr databases
# 7/8/2009 PW - revised fasta reader method to check amino acid characters
#             - added error checking flags to FasterReader & Protein methods
# 4/23/2010 PW - fixed findPeptide wrt I/L being indistinguishable
#              - improved adding decoy_string to accessions in reverseProtein
#              - added more support for Ref_Seq entries in NCBI databases
# 4/12/2011 PW - added a tryptic digest function
# 2011 PW - Switched to arrays to hold Gi to Taxon data to reduce memory use
# 2016 PW - rewritten for Python 3
# 2017 PW - removed any support for IPI databases and for cleaning accessions
#
import os
import sys
import gzip
import tarfile
import urllib.request
import socket
import sqlite3
import tkinter
from tkinter import filedialog

import fasta_lib

def get_folder(default_location, title_string=None):
    """Dialog box to browse to a folder.  Returns folder path.

    Usage: full_folder_name = get_folder(default_location, [title]),
        where "default_location" is a starting folder location,
        "title" is an optional message to list in the dialog box,
        and "full_folder_name" is the complete selected folder name.
    Written by Phil Wilmarth, 2008, 2016
    """
    # set up GUI elements
    root = tkinter.Tk()
    root.withdraw()
    try:
        root.tk.call('console', 'hide')
    except:
        pass
    
    # set default title string and location if not passed
    if title_string is None:   
        title_string = 'Select a folder with desired files/dirs'
    if not default_location:
        default_location = os.getcwd()
    
    # create dialog box for folder selection
    root.update()   # helps make sure dialog box goes away after selection
    full_folder_name = filedialog.askdirectory(parent=root, initialdir=default_location, 
                                               title=title_string, mustexist=True)    
    # return full folder name
    return full_folder_name   

def get_file(default_location, ext_list, title_string=None):
    """Dialog box to browse to a file.  Returns full file name.

    Usage: full_file_name = get_file(default_location, ext_list, [title]),
        where "default_location" is a starting folder location,
        ext_list is a list of (label, pattern) tuples,
        e.g. ext_list = [('Text files', '*.txt')],
        "title" is an optional message to list in the dialog box, and
        "full_file_name" is the complete name of the selected file.
    Written by Phil Wilmarth, OHSU, 2008, 2016.
    """
    # set up GUI elements
    root = tkinter.Tk()
    root.withdraw()
    try:
        root.tk.call('console', 'hide')
    except:
        pass
    
    # set default title string and ext list if not passed
    if title_string is None:   
        title_string = 'Select a single FILE'
    if not ext_list:
        ext_list =  [('All files', '*.*')]
    if not default_location:
        default_location = os.getcwd()
    
    # create dialog box for file selection
    root.update()   # helps make sure dialog box goes away after selection
    filename = filedialog.askopenfilename(parent=root, initialdir=default_location, 
                                          filetypes=ext_list, title=title_string)    
    # return full filename
    return filename       

def save_file(default_location, ext_list, default_file='', title_string=None):
    """Dialog box to save a file.  Returns full name of desired file.

    Usage: full_file_name = save_file(def_loc, ext_list, [def_file], [title]),
        where "def_loc" is a starting folder location,
        ext_list is a list of (label, pattern) tuples,
        e.g. ext_list = [('Text files', '*.txt')],
        "def_file" is an optional default filename,
        "title" is an optional message to list in dialog box, and
        "full_file_name" is the complete name of the desired file.
    Written by Phil Wilmarth, OHSU, 2009, 2016.
    """
    # set up GUI elements
    root = tkinter.Tk()
    root.withdraw()
    try:
        root.tk.call('console', 'hide')
    except:
        pass
    
    # set default title string if not passed
    if title_string is None:   
        title_string = 'Select a single FILE'
    if not ext_list:
        ext_list =  [('All files', '*.*')]
    if not default_location:
        default_location = os.getcwd()
        
    # create dialog box for file selection
    root.update()   # helps make sure dialog box goes away after selection
    filename = filedialog.asksaveasfilename(parent=root, initialdir=default_location, 
                                            initialfile=default_file, filetypes=ext_list, 
                                            title=title_string)
    # return full filename
    return filename   
    
def get_files(default_location, ext_list, title_string=None):
    """Dialog box to browse for files.  Returns a tuple of file names.

    Usage: file_name_list = get_files(default_location, ext_list, [title]),
        where "default_location" is a starting folder location,
        ext_list is a list of (label, pattern) tuples,
        e.g. ext_list = [('Text files', '*.txt')],
        "title" is an optional message to list in the dialog box, and
        "file_name_list" is a tuple of file name(s).
    Written by Phil Wilmarth, OHSU, 2010, 2016.
    """
    # set up GUI elements
    root = tkinter.Tk()
    root.withdraw()
    
    # set default title string if not passed
    if title_string is None:   
        title_string = 'Select one or more FILE(s)'
    if not ext_list:
        ext_list =  [('All files', '*.*')]
    if not default_location:
        default_location = os.getcwd()
        
    # create dialog box for file selection
    root.update()   # helps make sure dialog box goes away after selection
    filenames = filedialog.askopenfilenames(parent=root, initialdir=default_location, 
                                            filetypes=ext_list, multiple=True, 
                                            title=title_string)
    return filenames   

def get_string(title, prompt='Enter a string', initial=''):
    """Function to wrapper tkSimpleDialog.askstring function
    Written by Phil Wilmarth, OHSU, 2010.
    """
    from tkinter.simpledialog import askstring
    return askstring(title, prompt, initialvalue=initial)

def taxon_cmd_line_checker(argv):
    """Checks command line arguments for correctness.
    Returns dictionary of taxon, name pairs or empty dictionary.
    Written by Phil Wilmarth, OHSU, 2009.
    """
    tax_dict = {}
    if argv[0].endswith('.py'):
        argv = argv[1:]
    
    # need to have an even number of (integer, name) arguments
    try:
        pairs = [(argv[i], argv[i+1]) for i in range(0, len(argv), 2)]
        for (taxon, name) in pairs:
            for char in taxon:
                if not char.isdigit():
                    break
            tax_dict[abs(int(taxon))] = name
    
    # print usage information if error in format
    except:
        print('\n   ### invalid command line argument format ###\n')
        print('   arguments must be a series of "taxonomy number" "name" pairs')
        print('   where "taxonomy number" is integer and "name" is text string.')
        print('   example: prompt>python program.py 9606 human 4932 yeast\n')
    
    # return dictionary, it will be empty if any errors encountered
    return tax_dict

def string_cmd_line_checker(argv):
    """Checks command line arguments for correctness.
    Returns dictionary of string, name pairs or empty dictionary.

    Written by Phil Wilmarth, OHSU, 2009.
    """
    str_dict = {}
    if argv[0].endswith('.py'):
        argv = argv[1:]
    
    # need to have an even number of (string, name) arguments
    try:
        pairs = [(argv[i], argv[i+1]) for i in range(0, len(argv), 2)]
        for (string, name) in pairs:
            str_dict[string] = name
    
    # print usage information if error in format
    except:
        print('\n   ### invalid command line argument format ###\n')
        print('   arguments must be a series of "string" "name" pairs')
        print('   where "string" is a text pattern and "name" is text string.')
        print('   Note: enclose "string" in double quotes if it has whitespace.')
        print('   example: prompt>python program.py "Homo sapiens" human\n')
    
    # return dictionary, it will be empty if any errors encountered
    return str_dict


class Peptide:
    """Data structure for some basic peptide information
    """
    def __init__(self, sequence='', begin=0, end=0, mass=0, missed=0):
        self.seq = sequence
        self.beg = begin
        self.end = end
        self.mass = mass
        self.missed = missed
        return

class Protein:
    """Object to hold protein accession numbers, descriptions, and sequences.
    Methods:
        __init_:standard constructor, no parameters.
        readProtein: returns next protein from "fasta_reader"
        printProtein: prints sequence in FASTA format
        parseNCBI: cleans up nr entries
        parseUniProt: cleans up Sprot/Trembl entries
        parseCONT: cleans up Contaminant entries
        reverseProtein: reverses sequences and modifies accession/descriptions
        molwtProtein: computes average MW of sequence
        frequencyProtein: returns aa composition dictionary
        seqlenProtein: returns aa sequence length
        findPeptide: finds location of peptide in protein sequence
        coverage: calculates coverage and aa counts from peptide list        
    Written by Phil Wilmarth, OHSU, 2009.
    """
    def __init__(self):
        """Basic constructor, no parameters.
        """
        # bare bones __init__ function
        self.accession = 'blank'
        self.new_acc = 'blank'
        self.description = 'blank'
        self.new_desc = 'blank'
        self.sequence = ''
        self.sequence_padded = None
        self.sequence_masked = None
        self.pad_count = None
        self.length = 0
        self.peptides = []
        return

    def readProtein(self, fasta_reader):
        """Gets the next FASTA protein entry from FastaReader object.
        Usage: Boolean = object.readProtein(fasta_reader),
            where "object" is an instance of a Protein object and
            "fasta_reader" is an instance of a FastaReader object.
            Return value is "False" when EOF encountered.
        Written by Phil Wilmarth, OHSU, 2009.
        """
        status = fasta_reader.readNextProtein(self)
        self.new_acc = self.accession
        self.new_desc = self.description
        return status

    def printProtein(self, file_obj=None, length=80):
        """Prints FASTA protein entry to file (stdout is default).
        Usage: object.printProtein([file_obj=None, length=80]),
            where "object" is an instance of a Protein object, and
            "file_obj" is a file object (a value of None will print
            to standard out stream.  Optional "length" is number of
            characters per line for the protein sequence.
        Written by Phil Wilmarth, OHSU, 2009.
        """
        if file_obj == None:
            file_obj = sys.stdout
            
        # print new accession and new descriptor on first line
        if self.new_desc == '':
            print('>'+self.new_acc, file=file_obj)
        else:
            print('>'+self.new_acc, self.new_desc, file=file_obj)
        
        # initialize some things
        char_count = 0
        char_line = ''
        
        # build up sequence line with "length" characters per line
        for char in self.sequence:
            if char_count < length:  # do not have "width" chars yet
                char_line += char
                char_count += 1
            else:                   # line is "width" long so print and reset
                print(char_line, file=file_obj)
                char_line = char
                char_count = 1
        
        # print last sequence line (often less than "width" long) and return
        if len(char_line):
            print(char_line, file=file_obj)
        return

    def parseNCBI(self, REF_SEQ_ONLY=False):
        """Parses NCBI nr database accession numbers and descriptions
        Written by Phil Wilmarth, OHSU, 2009.
        """
        #===================================================
        # Might need to check this with new NCBI formats
        # maybe a regex would be more elegant if doable
        #===================================================
        # keep the gi number (or RefSeq) for the accession, get rid of others
        temp = self.accession.split('|')
        if REF_SEQ_ONLY:
            try:
                iacc = temp.index('ref')
            except ValueError:
                iacc = temp.index('gi')
        else:
            iacc = temp.index('gi')
        self.new_acc = '|'.join(temp[iacc:iacc+2])
        
        # keep the first description string if more than one
        new_desc = self.description.split(chr(1))[0]
        new_desc = new_desc.rstrip()
        
        # add a period to the end of description
        if not new_desc.endswith('.'):
            new_desc = new_desc + '.'
        self.new_desc = new_desc
        return

    def parseUniProt(self, KEEP_UNIPROT_ID=False):
        """Parses UniProt database accession numbers and descriptions.
        Written by Phil Wilmarth, OHSU, 2009.
        """
        if len(self.accession.split('|')) < 3:
            print('   parseUniProt WARNING: fewer than 3 accession elements')
        
        # if KEEP_UNIPROT_ID is set to True, keep the more human readible part of the accession numbers
        if KEEP_UNIPROT_ID:
            self.new_acc = self.accession.split('|')[-1]
            acc_num = self.accession.split('|')[-2]
        else:
            self.new_acc = self.accession.split('|')[-2]
            acc_num = self.accession.split('|')[-1]
        
        # keep the desciption text and get rid of the trailing stuff
        new_desc = self.description.split(chr(1))[0]
        new_desc = new_desc.split('OS=')[0]
        new_desc = new_desc.rstrip()
        
        # add other accession text to end of description
        if new_desc.endswith('.'):
            new_desc = new_desc[:-1]
        self.new_desc = '%s (%s).' % (new_desc, acc_num)
        return

    def parseCONT(self):
        """Parses contaminant accession numbers and descriptions
        Written by Phil Wilmarth, OHSU, 2009.
        """
        # keep the contaminant tag for the accession, get rid of others
        old_acc = ''
        temp = self.accession.split('|')
        if len(temp) > 1:
            self.new_acc = temp[0]
            old_acc = '|'.join(temp[1:])
        
        # keep the first description string if more than one
        new_desc = self.description.split(chr(1))[0]
        new_desc = new_desc.strip()
        
        # add other accessions to description (if any)
        if new_desc.endswith('.'):
            new_desc = new_desc[:-1]
        if old_acc:
            self.new_desc = '%s (%s).' % (new_desc, old_acc)
        else:
            self.new_desc = new_desc + '.'
        return

    def reverseProtein(self, decoy_string):
        """Reverses protein sequence and returns new Protein object.
        Usage: rev_prot = object.reverseProtein(decoy_string),
            where "object" is a Protein object, "decoy_string" is the
            unique identifier text to add to the beginning of the 
            protein accesion number, and "rev_prot" is new Protein object.
        Written by Phil Wilmarth, OHSU, 2009.
        """
        # make sure decoy_string ends with an undescore
        if not decoy_string.endswith('_'):
            decoy_string = decoy_string + '_'
        
        # create a new Protein instance
        rev_prot = Protein() 
        
        # prefix the decoy_sting to all important parts of accession
        temp = self.accession.split('|')
        if self.accession.startswith('gi|'):    # modify nr DB accessions
            temp = [temp[i] for i in range(1, len(temp), 2)]
            temp = temp[0:1] # keep the gi number
        elif self.accession.startswith('sp|') or \
             self.accession.startswith('tr|'):  # modifiy Sprot or Trembl accs
            temp = temp[1:]
            temp = temp[0:1] # keep the accession number
        elif self.accession.startswith('CONT_'):
            temp = temp[0:1] # keep the contaminant number
        else:
            pass
        #
        if len(temp) > 1: # make sure temp is a list not a string
            rev_prot.accession = '|'.join([decoy_string+x for x in temp])
        else:
            rev_prot.accession = decoy_string + temp[0]
        rev_prot.new_acc = rev_prot.accession

        # change the desciptions, too.
        rev_prot.description = 'REVERSED.'
        rev_prot.new_desc = 'REVERSED.'
        
        # now reverse the sequence
        rev_prot.sequence = self.sequence[::-1]
        return rev_prot

    def molwtProtein(self, show_errs=True):
        """Returns protein molecular weight as the sum of average aa masses.
        If "show_errs" flag set, invalid amino acid characters are reported.

        Written by Phil Wilmarth, OHSU, 2009.
        """        
        # start with water and H+ masses, then add aa masses
        bad_char = {}
        self.setMasses()
        molwt = 18.01 + 1.007825
        for i in self.sequence:
            try:
                molwt += self.ave_masses[i]
            except:     # keep track of bad characters
                bad_char[i] = True
        #
        bad_char = sorted(list(bad_char.keys()))
        if len(bad_char) > 0 and show_errs:     # report bad chars if desired
            print('   WARNING: unknown symbol(s) (%s) in %s:\n%s' %
                  (''.join(bad_char), self.accession, self.sequence))
        return molwt

    def frequencyProtein(self, show_errs=True):
        """Returns aa frequency distrubution as a dictionary.
        If "show_errs" flag set, invalid amino acid characters are reported.        
        Written by Phil Wilmarth, OHSU, 2009.
        """
        freq = {'X':0, 'G':0, 'A':0, 'S':0,
                'P':0, 'V':0, 'T':0, 'C':0,
                'L':0, 'I':0, 'J':0, 'N':0,
                'O':0, 'B':0, 'D':0, 'Q':0,
                'K':0, 'Z':0, 'E':0, 'M':0,
                'H':0, 'F':0, 'R':0, 'Y':0,
                'W':0, 'U':0, '*':0, '-':0}
        
        # count the amino acids for all residues in sequence
        bad_char = {}
        for i in self.sequence:
            try:
                freq[i] += 1
            except:     # keep track of bad characters
                bad_char[i] = True
        #
        bad_char = list(bad_char.keys())
        bad_char.sort()
        if len(bad_char) > 0 and show_errs: # report any bad chars, if desired
            print('   WARNING: unknown symbol(s) (%s) in %s:\n%s' %
                  (''.join(bad_char), self.accession, self.sequence))
        return freq

    def seqlenProtein(self):
        """Calculates protein sequence length.
        Written by Phil Wilmarth, OHSU, 2009.
        """
        return (len(self.sequence))

    def findPeptide(self, peptide, pad_count=1):
        """Calculates location of all 'peptide' matches in 'self.sequence.'

        Written by Phil Wilmarth, OHSU, 2009, 2016.
        """
        import re
        matches = []
        
        # get rid of bounding residues, if any
        try:
            peptide = peptide.split('.')[1]
        except IndexError:
            pass
        
        # remove any modification symbols and mask I/L:
        # '*', '#', '@', '^', '~', '$', '%', '!', '+', 'n', 'c', '[', ']' (Current Comet along with old style nt, ct)
        splitter = re.compile(r'[*#@^~$%!+nc\[\]]')
        base_pep = ''.join(splitter.split(peptide))
        base_pep_masked = re.sub(r'[IL]', 'j', base_pep)
 
        # fix the protein sequence for peptide lookups (pad and mask I/L). Save the results to improve performance
        if (not self.sequence_masked) or (pad_count != self.pad_count):
            self.sequence_padded = ('-' * pad_count) + self.sequence + ('-' * pad_count) # add bounding symbols
            self.sequence_masked = re.sub(r'[IL]', 'j', self.sequence_padded)
            self.pad_count = pad_count
        
        # find all matches of base_pep_masked to protein sequence (padded and masked)
        for match in re.finditer(base_pep_masked, self.sequence_masked):
            start, end = match.span()
            start_prot, end_prot = start - self.pad_count + 1, end - self.pad_count

            # add bounding AAs, periods, and put back modification special chars
            pre = self.sequence_padded[start-self.pad_count:start]
            post = self.sequence_padded[end:end+self.pad_count]
            full_seq = pre + '.' + peptide + '.' + post
            matches.append((start_prot, end_prot, full_seq))
        
        # return the match list (empty list if no matches)
        return matches

    def calcCoverage(self, peptide_list):
        """Calculates % coverage and aa frequency map of matched peptides.
        "peptide_list" is list of sequences with optional counts (as tuples).
        Written by Phil Wilmarth, OHSU, 2009.
        """
        freq_dict = {}
        try:    # see if peptide_list is a list of tuples or not
            for peptide, count in peptide_list:
                for (beg, end, seq) in self.findPeptide(peptide):
                    for key in [str(i) for i in range(beg, end+1)]:
                        if freq_dict.get(key, False):
                            freq_dict[key] = freq_dict[key] + count
                        else:
                            freq_dict[key] = count
        except ValueError:
            for peptide in peptide_list:
                for (beg, end, seq) in self.findPeptide(peptide):
                    for key in [str(i) for i in range(beg, end+1)]:
                        if freq_dict.get(key, False):
                            freq_dict[key] = freq_dict[key] + 1
                        else:
                            freq_dict[key] = 1
        #
        coverage = 100.0*float(len(freq_dict))/float(len(self.sequence))
        coverage_map = []
        for i, aa in enumerate(self.sequence):
            coverage_map.append((str(i+1), aa, freq_dict.get(str(i+1), 0)))
        return (coverage, coverage_map)
    
    def enzymaticDigest(self, enzyme_regex=None, low=500.0, high=5000.0, length=7, missed=2, mass='mono'):
        """Performs a tryptic digest of a protein sequence. This does not
        do any modifications to residues except for reduction/alkylation of
        cys residues (C+57). Mass filters should be relaxed.
        
        Returns a list of digested peptides.
        enzyme_regex is a compiled re object for the enzyme cleavage
            (if enzyme_regex not defined, do tryptic digest by default)
        low, high - mass limits for peptides.
        length - minimum amino acid length
        missed - maximum number of missed cleavages.
        mass - 'ave' average or 'mono' monoisotopic masses.

        written by Phil Wilmarth, OHSU, 2011, 2016.
        """
        """Regular expression digestion table:
        trypsin from: http://stackoverflow.com/questions/18364380/python-3-cut-peptide-regular-expression

        regex = re.compile(r".")                        # no enzyme
        regex = re.compile(r".(?:(?<![KR](?!P)).)*")    # trypsin strict
        regex = re.compile(r".(?:(?<![KR]).)*")         # trypsin with cleavage at P
        regex = re.compile(r".(?:(?<![K](?!P)).)*")     # Lys-C strict
        regex = re.compile(r".(?:(?<![K]).)*")          # Lys-C with cleavage at P
        regex = re.compile(r".(?:(?![K]).)*")           # Lys-N
        regex = re.compile(r".(?:(?<![R](?!P)).)*")     # Arg-C strict
        regex = re.compile(r".(?:(?![D]).)*")           # Asp-N
        regex = re.compile(r".(?:(?<![M]).)*")          # CnBr
        regex = re.compile(r".(?:(?<![DE](?!P)).)*")    # Glu-C
        regex = re.compile(r".(?:(?<![FL](?!P)).)*")    # PepsinA
        regex = re.compile(r".(?:(?<![FWYL](?!P)).)*")  # chymotrypsin
        """   
        import copy
        import re

        # skip if there is no sequence to digest
        if len(self.sequence) == 0:
            return []

        # tryptic digestion is the default
        if not enzyme_regex:
            enzyme_regex = re.compile(r".(?:(?<![KR](?!P)).)*")
        
        # set up masses, default is alkylated cysteine. No mechanism for other modifications yet.
        self.setMasses()
        if mass == 'ave':
            masses = copy.deepcopy(self.ave_masses)
            masses['C'] = 160.197
        elif mass == 'mono':
            masses = copy.deepcopy(self.mono_masses)
            masses['C'] = 160.03065
        else:
            print('...WARNING: masses must be "ave" or "mono"')

        # digest the sequence
        digest_matches = [x for x in enzyme_regex.finditer(self.sequence)] # list of re match objects

        # get info from match objects into Peptide object attributes
        digest = [Peptide(mass=masses['water']) for x in digest_matches]
        for i, match in enumerate(digest_matches):
            digest[i].seq = match.group()
            digest[i].beg, digest[i].end = match.span()
            digest[i].beg += 1
            for aa in match.group():
                digest[i].mass += masses[aa]
            
        # test peptides and missed cleavage peptides for mass ranges and min length
        valid_digest = []
        for i in range(len(digest)):
            
            # check if peptide is within the mass range and meets min length
            if (low <= digest[i].mass <= high) and (len(digest[i].seq) >= length):
                valid_digest.append(digest[i])
                
            # create and check missed cleavages
            for j in range(1, missed+1):
                if (i+j) > len(digest)-1:
                    continue
                temp = Peptide(begin=100000)    # a peptide object for missed cleavages
    
                # calculate running sums for each number of missed cleavages
                for k in range(j+1):
                    if (i+k) > len(digest)-1:
                        continue
                    temp.seq += digest[i+k].seq
                    temp.beg = min(temp.beg, digest[i+k].beg)
                    temp.end = max(temp.end, digest[i+k].end)
                    temp.mass += (digest[i+k].mass - masses['water'])
                temp.mass += masses['water']
                temp.missed = k
                
                # check missed cleavage peptide for valid mass range and length
                if (low <= temp.mass <= high) and (len(temp.seq) >= length):
                    valid_digest.append(temp)
        
        # return the list of digested peptides
        self.peptides = valid_digest
        return valid_digest

    def setMasses(self):
        """Set average and monoisotopic mass dictionaries.
        """
        self.ave_masses = {'X':  0.0000, 'G': 57.0513, 'A': 71.0779, 'S': 87.0773, 'P': 97.1152,
                           'V': 99.1311, 'T':101.1039, 'C':103.1429, 'L':113.1576, 'I':113.1576,
                           'J':113.1576, 'N':114.1026, 'O':114.1472, 'B':114.5950, 'D':115.0874,
                           'Q':128.1292, 'K':128.1723, 'Z':128.6216, 'E':129.1140, 'M':131.1961,
                           'H':137.1393, 'F':147.1739, 'R':156.1857, 'Y':163.1733, 'W':186.2099,
                           'U':150.0379, '*': 0.00000, '-': 0.00000, 'water':18.02}
        self.mono_masses = {'X':  0.000000, 'G': 57.021464, 'A': 71.037114, 'S': 87.032028, 'P':97.052764,
                            'V': 99.068414, 'T':101.047679, 'C':103.009185, 'L':113.084064, 'I':113.084064,
                            'J':113.084064, 'N':114.042927, 'O':114.147200, 'B':114.595000, 'D':115.026943,
                            'Q':128.058578, 'K':128.094963, 'Z':128.621600, 'E':129.042593, 'M':131.040485,
                            'H':137.058912, 'F':147.068414, 'R':156.101111, 'Y':163.063320, 'W':186.079313,
                            'U':150.953630, '*':  0.000000, '-':  0.000000, 'water':18.01057}
        return

    # end class

class FastaReader:
    """Reads FASTA entries from a file-like object.
    methods:
    __init__: basic constructor, no parameters.    
    readProtein: reads one FASTA entry from a file object (text or zipped)
        arguments are "next_protein" and "file_obj"
        returns True (next protein) or False (EOF or not FASTA).
    written by Phil Wilmarth, OHSU, 2009.
    """

    def __init__(self, fasta_file):
        """Basic constructor function.  No parameters
        self._last_line retains the previous '>' line and
        self._valid is a dictionary of valid protein FASTA chars.
        """
        # attribute to save last line from previous read
        self._last_line = 'start value'
        self._file_obj = None
        self._fasta_file = fasta_file
        
        # list of valid amino acid characters
        self._valid = {'X':True, 'G':True, 'A':True, 'S':True, 'P':True,
                       'V':True, 'T':True, 'C':True, 'L':True, 'I':True,
                       'J':True, 'N':True, 'O':True, 'B':True, 'D':True,
                       'Q':True, 'K':True, 'Z':True, 'E':True, 'M':True,
                       'H':True, 'F':True, 'R':True, 'Y':True, 'W':True,
                       'U':True, '*':True, '-':True }

##        if not os.path.exists(fasta_file):
##            ext_list = [('FASTA files', '*.fasta'), 
##                        ('Zipped FASTA files', '*.gz'), 
##                        ('All files', '*.*')]
##            fasta_file = get_file(os.getcwd(), ext_list, title_string="Select FASTA file")
        
        # get file object and save as attribute
        try:
            if fasta_file.endswith('.gz'):
                self._file_obj = gzip.open(fasta_file, 'rt')
            else :
                self._file_obj = open(fasta_file, 'Ur')
        except IOError:
            print('   WARNING:', fasta_file, 'could not be opened!')
            raise
        return

    def readNextProtein(self, next_protein, check_for_errs=False):
        """Loads one FASTA protein text entry into a Protein object.
        Returns True (protein entry found) or False (end of file).
        If "check_for_errs" flag is set, amino acid chars are checked.
        Written by Phil Wilmarth, OHSU, 2009.
        """
        # at first call, start reading lines
        if self._last_line == 'start value':
            self._last_line = self._file_obj.readline()
            if not self._last_line:
                self._file_obj.close()
                return(False)
            self._last_line = self._last_line.strip()
        
        # get next protein's info from _last_line
        if self._last_line.startswith('>'):
            next_protein.accession = self._last_line.split()[0][1:]
            next_protein.new_acc = next_protein.accession
            start = len(next_protein.accession)+2
            next_protein.description = self._last_line[start:]
            next_protein.new_desc = next_protein.description
        
        # return if empty line (EOF) or non-description line
        else:
            self._file_obj.close()
            return(False)                    
        
        # reset variables and read in next entry
        next_protein.sequence = ""
        line = self._last_line
        self._last_line = ""
        bad_char = {}
        while line:
            line = self._file_obj.readline()
            if not line:
                break
            else:
                testline = line.strip()
            if testline == '':
                continue
            
            # stop reading at next descriptor line (and save line)
            if line.startswith('>'):
                self._last_line = line.strip()
                
                # report bad characters if conditions were met
                bad_char = sorted(bad_char.keys())
                if len(bad_char) > 0 and check_for_errs:
                    print('   WARNING: unknown symbol(s) (%s) in %s' %
                          (''.join(bad_char), next_protein.accession))
                break
            
            # add next sequence line to protein's sequence
            else:
                line = line.rstrip()
                line = line.upper()
                if check_for_errs: # checking chars slows down the program
                    for char in line:
                        if self._valid.get(char, False):
                            next_protein.sequence += char
                        else:
                            bad_char[char] = True                
                else: # blindly adding the line is faster...
                    next_protein.sequence += line
                    
        # return (protein info retained in next_protein)
        return True
        
    # end class

    
def get_uniprot_version():
    """Gets UniProt version numbers from online release notes.
    Written by Phil Wilmarth, OHSU, 2009.
    """
    # set up to read the online UniProt release notes file
    print('...getting database version numbers...')
    versions = {'uniprot':'XX.X', 'sprot':'XX.X', 'trembl':'XX.X'}
    address = 'ftp://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/complete/reldate.txt'
    reldate = urllib.request.urlopen(address)
    charset = 'utf-8' # should try an get this dynamically (I was getting None with header calls)
    
    # read release notes file and get UniProt, etc. version numbers
    for line in reldate:
        line = line.decode(charset)
        if 'UniProt Knowledgebase' in line:
            versions['uniprot'] = line.split()[3].replace('_', '.')
        elif 'Swiss-Prot' in line:
            versions['sprot'] = line.split()[2].replace('_', '.')
        elif 'TrEMBL' in line:
            versions['trembl'] = line.split()[2].replace('_', '.')
    return(versions)

def download_uniprot(db, folder, versions):
    """Downloads (if necessary) UniProt FASTA, taxon files to "folder".
    Written by Phil Wilmarth, OHSU, 2009.
    """
    # check if files are already downloaded, if not fetch them from uniprot site
    socket.setdefaulttimeout(120.)
    print('...downloading databases and taxonomy files...')
    db_name = 'uniprot_%s_%s.fasta.gz' % (db, versions[db],)
    base_address = 'ftp://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/complete/'
    db_address = 'uniprot_%s.fasta.gz' % (db,)
    db_address = base_address + db_address
    files_addresses = [ (os.path.join(folder, db_name), db_address),
                         (os.path.join(folder, 'taxdump.tar.gz'),
                         'ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz'),
                         (os.path.join(folder, 'speclist.txt'),
                         'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/complete/docs/speclist.txt')]
    files_addresses.reverse()
    for (file_name, address) in files_addresses:
        if os.path.exists(file_name):
            pass
        else:
            print('...downloading', file_name)
            x = reporter()
            try:
                urllib.request.urlretrieve(address, file_name, reporthook=x.report)
            except:
                print('...WARNING: download may have hung at EOF or other error')
            finally:
                urllib.request.urlcleanup()
    return

def download_ncbi(nr_folder):
    """Downloads (if necessary) the ncbi FASTA and taxon files to "nr_folder".
    Written by Phil Wilmarth, OHSU, 2009.
    """

    # check if files are already downloaded, if not fetch them from ncbi site
    nr_name = os.path.split(nr_folder)[1] + '.gz'
    if os.path.exists(os.path.join(nr_folder, 'nr.gz')):
        os.rename(os.path.join(nr_folder, 'nr.gz'), \
                  os.path.join(nr_folder, nr_name))
    files_addresses = [ (os.path.join(nr_folder, nr_name),
                         'ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nr.gz'),
                        (os.path.join(nr_folder, 'prot.accession2taxid.gz'),
                         'ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz'),
                        (os.path.join(nr_folder, 'taxdump.tar.gz'),
                         'ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz') ]
    files_addresses.reverse()
    socket.setdefaulttimeout(120.)
    for (file_name, address) in files_addresses:
        if os.path.exists(file_name):
            pass
        else:
            print('...downloading', file_name)
            x = reporter()
            try:
                urllib.request.urlretrieve(address, file_name, reporthook=x.report)
            except:
                print('...WARNING: download may have hung at EOF or other error')
            finally:
                urllib.request.urlcleanup()
    return


class reporter():
    """Prints download progress to console.
    """
    def __init__(self):
        self.packets = 0
        self.size = 0
        self.buff = 8192
        
    def report(self, packets, buff, size):
        # this monitors the download progress
        if packets % 512 == 0:
            sub_total = packets * buff
            print('......%s of %s bytes (%.2f%%)' %
                  ("{0:,d}".format(sub_total), "{0:,d}".format(size), float(100*sub_total)/size))
        return

    # end class


def expand_species(folder, db, taxon_dict, min_sequence_count, min_seq_per_species,
                   REF_SEQ_ONLY=False):
    """Expands any taxon nodes numbers into all member taxon numbers.
    Written by Phil Wilmarth, OHSU, 2009.
    """
    VERBOSE = False
    
    # open taxonomy nodes file
    print('...making taxonomy nodes dictionary...')
    archive_name = os.path.join(folder, 'taxdump.tar.gz')
    archive = tarfile.open(archive_name)
    nodes = archive.extractfile('nodes.dmp')

    # read nodes file and save child to parent taxon mappings
    child_to_parent = {}
    while True:
        line = nodes.readline()
        line = line.decode('utf-8')
        line = line.rstrip()
        if not line:
            break
        else:
            line = line.rstrip()
        item = line.split('\t|\t')
        child_to_parent[int(item[0])] = int(item[1])    
    nodes.close()
    
    # open the fasta_analysis.txt file
    species_counts = {}
    if db == 'nr':
        analysis_file = os.path.join(folder, 'nr_fasta_analyze.txt')
        if REF_SEQ_ONLY:
            index = 4
        else:
            index = 3
    elif db == 'sprot':
        analysis_file = os.path.join(folder, 'sprot_fasta_analyze.txt')
        index = 4
    elif db == 'trembl':
        analysis_file = os.path.join(folder, 'trembl_fasta_analyze.txt')
        index = 4
    else:
        analysis_file = os.path.join(folder, 'uniprot_fasta_analyze.txt')
        index = 6
    fasta_analyze = open(analysis_file, 'r')
    
    # save species taxons and sequence counts
    line = fasta_analyze.readline().rstrip()    # skip header line
    while True:
        line = fasta_analyze.readline()
        if not line:
            break
        else:
            line = line.rstrip()    
        temp = line.split('\t')
        taxon = temp[1]
        count = temp[index]
        try:
            taxon = int(taxon)
            count = int(count)
        except:
            continue
        species_counts[taxon] = count
    fasta_analyze.close()
    
    # see if we have any group taxon numbers
    group_expand = {}
    for alltax in [x for x in species_counts.keys() if species_counts[x] >= min_seq_per_species]:
        tree = []
        tree.append(alltax)
        parent = alltax
        while parent != 1:
            try:
                parent = child_to_parent[parent]
            except KeyError:
                if VERBOSE:
                    print('...WARNING: lookup of taxon %s failed' % (parent,))
                break
            tree.append(parent)
        
        # see if any of our taxon dictionary numbers are in lineage
        for tax in list(taxon_dict.keys()):
            if (tax in tree) and (alltax != tax):
                    taxon_dict[alltax] = taxon_dict[tax]
                    try:
                        group_expand[tax].append(alltax)
                    except KeyError:
                        group_expand[tax] = [alltax]
    
    # report any expanded taxononmy numbers
    for tax in taxon_dict.keys():
        if len(group_expand.get(tax, [])) >= 1:
            print('...taxon %s (%s direct seqs) was a node with children:'
                  % (tax, "{0:,d}".format(species_counts.get(tax, 0))))
            for child in group_expand[tax]:
                print('......taxon %s has %s sequences' %
                      (child, "{0:,d}".format(species_counts[child])))                
    
    # check taxon dictionary and remove taxons with too few sequences
    for tax in list(taxon_dict.keys()):
        if species_counts.get(tax, 0) <= min_sequence_count:
            print('...WARNING: taxon number %s had too few proteins' % (tax,))
            del taxon_dict[tax]
    return

class AccToTaxon(): # rename this
    """Object to map NCBI accessions to taxonomy numbers.
    Methods:
        __init__: loads or arrays from taxonomy files (or reloads)
        get(acc, default): return taxon number of "acc" or "default"
    Written by Phil Wilmarth, OHSU, 2009,2017.
    """
    def __init__(self, nr_folder):
        # mostly placeholders
        self.conn = None    # connection to sqlite3 database
        self.which = 'sqlite3' # no longer needed flag
        return

    def create_or_load(self, nr_folder):
        """Creates or reloads an AccToTaxon object.
        """
        # build filename
        acc_to_tax_db = os.path.join(nr_folder, 'acc_to_tax_DB.sq3')
        
        # check for SQLITE3 database file existance
        if os.path.exists(acc_to_tax_db):
            print('...conecting to SQLite database...')
            self.conn = sqlite3.connect(acc_to_tax_db)

        # if DB file not found, create it
        else:
            print('...making acc_to_taxon mapping database...')
            
            # read the gi number to taxonomy number data into an SQLite DB
            self.conn = sqlite3.connect(acc_to_tax_db)    # create and connect to DB
            c = self.conn.cursor()
            c.execute('''CREATE TABLE acc_to_tax(
                            acc text PRIMARY KEY NOT NULL,
                            tax integer NOT NULL)''')
            line_count = 0
            fname = os.path.join(nr_folder, 'prot.accession2taxid.gz')
            for line in gzip.open(fname, mode='rt'):
                line = line.rstrip()
                if 'accession.version' in line:     # skip header line
                    continue
                line_count += 1
                if(line_count % 10000) == 0: # update databae periodically
                    self.conn.commit()
                if (line_count % 1000000) == 0:
                    print('......(%s acc_to_taxon lines read)' % ("{0:,d}".format(line_count),))
                    
                # parse the accession.version and taxonomy fields
                acc = line.split('\t')[0]
                tax = int(line.split('\t')[2])
                c.execute('INSERT INTO acc_to_tax (acc, tax) VALUES (?, ?)', (acc, tax))
            self.conn.commit()  # make sure databae is fully updated
            c.close()
        return
    
    def get(self, acc, default):
        """Lookup of taxon number given an NCBI accession.
        Switched to SQLite database to avoid out-of-memory errors
            Phil W., 2013-7.
        """
        # query database for taxon number associated with accession        
        c = self.conn.execute('SELECT tax FROM acc_to_tax WHERE acc=?', (acc,))
        try:
            tax = c.fetchone()[0]
        except TypeError:   # accession not in database
            tax = default
        return tax

    def close(self):
        """Closes SQLite database connection.
        written by Phil Wilmarth, OHSU, 2013
        """
        try:
            self.conn.close()
        except:
            pass
        return
    
    # end class

def make_taxon_to_sci_name(folder):
    """Makes the taxon_to_sci_name dictionary.
    Written by Phil Wilmarth, OHSU, 2009.
    """
    print('...making taxon_to_sci_name dictionary...')
    archive_name = os.path.join(folder, 'taxdump.tar.gz')
    archive = tarfile.open(archive_name)
    names = archive.extractfile('names.dmp')
    taxon_to_name = {}
    
    # read file and save names from 'scientific name' lines
    line = names.readline()
    line = line.decode('utf-8')
    line = line.rstrip()
    while line:
        item = line.split('\t')
        if item[6] == 'scientific name':
            taxon_to_name[int(item[0])] = item[2]
        line = names.readline()
        line = line.decode('utf-8')
        line = line.rstrip()
    names.close()
    
    # there may be some gi numbers that have taxon id of zero
    taxon_to_name[0] = 'Zero_taxon_number'
    return taxon_to_name

def make_uniprot_to_taxon(folder):
    """Makes sci_to_taxon and id_to_taxon dictionaries from "speclist.txt".
    Written by Phil Wilmarth, OHSU, 2009.
    """
    print('...making scientific names and IDs to taxon dictionaries...')
    speclist = open(os.path.join(folder, 'speclist.txt'), 'r')
    sci_to_taxon = {}
    id_to_taxon = {}
    
    # read file and save names from 'scientific name' lines
    start = False
    while True:
        line = speclist.readline()
        if not line: break
        if 'Real organism codes' in line:
            start = True
        if '"Virtual" codes' in line:
            start = False
        if start:
            line = line.rstrip()
            line = line.lstrip()
            item = line.split('N=')
            if len(item) > 1:
                name = item[1]
                if name == 'Official (scientific) name':
                    continue
                bit = item[0].split()
                spec_id = '_' + bit[0]
                taxon = bit[2][:-1]
                sci_to_taxon[name] = int(taxon)
                id_to_taxon[spec_id] = int(taxon)
    
    # close file and return dictionaries
    speclist.close()
    return sci_to_taxon, id_to_taxon

def make_all_names_to_taxon(folder):
    """Makes the all_names_to_taxon dictionary.
    Written by Phil Wilmarth, OHSU, 2009.
    """
    print('...making all_names_to_taxon dictionary...')
    archive_name = os.path.join(folder, 'taxdump.tar.gz')
    archive = tarfile.open(archive_name)
    names = archive.extractfile('names.dmp')
    all_names_to_taxon = {}
    
    # read file and save taxonomy numbers for all names
    while True:
##        line = names.readline().rstrip()
        line = names.readline()
        line = line.decode('utf-8')
        line = line.rstrip()
        if not line: break          
        item = line.split('\t')
        name = item[2].replace('"','')
        name = name.lstrip()
        all_names_to_taxon[name] = int(item[0])
    
    # close archive and return dictionary
    names.close()
    return all_names_to_taxon

def uniprot_species_frequency(database_name):
    """Compiles species frequency info from Sprot or Trembl databases.
    Written by Phil Wilmarth, OHSU, 2009.
    """
    # read all of the protein descriptions and parse out species names
    if 'sprot' in database_name:
        print('...scanning sprot database for species names...')
    else:
        print('...scanning trembl database for species names (slow)...')
    name_freq = {}
    name_to_spec_id = {}
    prot_count = 0
    f = gzip.open(database_name, 'rt')
    while True:
        line = f.readline()
        if not line: break
        
        # get species name, id; save in dictionary; make frequency totals
        if line.startswith('>'):
            prot_count += 1
            if (prot_count % 500000) == 0:
                print('......(%s proteins read...)' % ("{0:,d}".format(prot_count),))
            (spec_id, name) = uniprot_parse_line(line)
            name_to_spec_id[name] = spec_id
            fasta_lib.add_or_increment(name, name_freq)            
    return name_freq, name_to_spec_id, prot_count

def uniprot_parse_line(line):
    """Parses UniProt description lines for species IDs and names.
    Written by Phil Wilmarth, OHSU, 2009.
    """
    name = 'Unknown species'
    spec_id = '_Z1234'
    
    # get the species ID
    spec_id = line.split()[0]
    spec_id = spec_id.split('|')[-1]
    spec_id = spec_id.split('_')[1]
    spec_id = '_' + spec_id
    
    # get the species name if present
    if len(line.split('OS=')) == 1:
        pass
    else:
        name = line.split('OS=')[1]
        if len(name.split('=')) == 1:
            name = name.rstrip()
        else:
            name = name.split('=')[0][:-2]
            name = name.rstrip()
    return spec_id, name           

def add_or_increment(species_name, species_dict):
    """Increments species count (if existing) or adds species to dictionary.
    Written by Phil Wilmarth, OHSU, 2009.
    """
    # add new species name or increment counter if name already exists
    try:
        species_dict[species_name] += 1
    except:
        species_dict[species_name] = 1
    return

def sort_species(species_dict):
    """Sorts a dictionary of species frequencies by genus name then frequency.
    Written by Phil Wilmarth, OHSU, 2009.
    """
    # first sort by Genus name using the "sort" attribute
    sorted_list = list(species_dict.items())
    sorted_list.sort()
    return sorted_list
##    
##    # sort the Specific names by inverse protein sequence count
##    # we need to define the Genus subsections of the database
##    groups = []
##    current = sorted_list[0][0]
##    start = 0
##    end = 0
##    
##    # find out where each Genus starts and stops in the list
##    for i, species in enumerate(sorted_list):
##        end += 1
##        if species[0] != current:
##            groups.append((start, end-1))
##            start = end-1
##            current = species[0]
##    groups.append((start, end))
##    
##    # now we can sort by sequence count ("second_sort" attribute) and then reverse for descending order
##    for (start, end) in groups:
##        temp = [(y, x) for (x, y) in sorted_list[start:end]]
##        temp.sort()
##        temp.reverse()
##        sorted_list[start:end] = [(y, x) for (x, y) in temp]
##    return sorted_list

def virus_test(name):
    """Tests whether "virus" is in a species name or not.
    Written by Phil Wilmarth, OHSU, 2009.
    """
    return ('virus' in name)

def save_species_info_nr(folder, name_freq, name2tax, ref2freq=None, minimum=0):
    """Writes taxonomy and species name info to file.
    Written by Phil Wilmarth, OHSU, 2009.
    """
    # sort the species names and write to file
    print('...writing species analysis results file...')
    sorted_list = fasta_lib.sort_species(name_freq)
    fout_name = os.path.join(folder, 'nr_fasta_analyze.txt')
    fout = open(fout_name, 'w')
    sformat = '%s\t%s\t%s\t%s\t%s\t%s\t%s'
    print(sformat % ('=SUBTOTAL(109,A2:A65000)', 'Taxon_ID', 'Species_Name', 'Sequence_Count',
                     'RefSeq_Count', 'Word_Count', 'Is_virus'), file=fout)
    dict_list = [name2tax, None, None, None]
    for name, count in sorted_list:
        if int(count) >= minimum:
            taxon = fasta_lib.get_taxon_from_name('nr', name, dict_list)
            ref_count = ref2freq.get(taxon, '')
            print(sformat % (str(1), taxon, name, count, ref_count, len(name.split()),
                             fasta_lib.virus_test(name)), file=fout)
    fout.close()
    return

def save_species_info(db, folder, name_freq, name2tax, sci2tax=None,
                      id2tax=None, name2id=None, minimum=0):
    """Writes taxonomy and species name info to file.
    Written by Phil Wilmarth, OHSU, 2009.
    """
    # sort the species names and write to file
    print('...writing species analysis results file...')
    sorted_list = fasta_lib.sort_species(name_freq)
    fout_name = os.path.join(folder, db+'_fasta_analyze.txt')
    fout = open(fout_name, 'w')
    sformat = '%s\t%s\t%s\t%s\t%s\t%s\t%s'
    print(sformat % ('=SUBTOTAL(109,A2:A65000)', 'Taxon_ID', 'Species_ID', 'Species_Name',
                     'Sequence_Count', 'Word_Count', 'Is_virus'), file=fout)
    dict_list = [name2tax, sci2tax, id2tax, name2id]
    for name, count in sorted_list:
        if int(count) >= minimum:
            taxon = fasta_lib.get_taxon_from_name(db, name, dict_list)
            #
            try:
                ID = name2id[name]
            except:
                ID = ' '
            print(sformat % (str(1), taxon, ID, name, count, len(name.split()),
                             fasta_lib.virus_test(name)), file=fout)                
    fout.close()
    return

def combine_analysis_files(folder):
    """Program to add Sprot columns to Trembl analysis file
    Written by Phil Wilmarth, OHSU, 2009.
    """
    trembl = open(os.path.join(folder, 'trembl_fasta_analyze.txt'), 'r')
    sprot = open(os.path.join(folder, 'sprot_fasta_analyze.txt'), 'r')
    
    header = trembl.readline().rstrip()
    new_header = header.split('\t')
    new_header[4] = 'Trembl_Count'
    new_header.insert(4, 'Sprot_Count')
    new_header.insert(6, 'Total_Count')
    header = '\t'.join(new_header)
    
    trembl_dict = {}
    i = 0
    while True:
        trline = trembl.readline()
        i += 1
        if not trline:
            break
        else:
            trline = trline.rstrip()
        trline = trline.split('\t')
        trline.insert(4, '')
        trline.insert(6, trline[5])
        tax = trline[1]
        if int(trline[5]) >= 5:
##            trembl_dict[tax] = [trline[3].split()[0], trline]
            trembl_dict[tax] = [trline[3], trline]
    sprot.readline() # skip header
    missing = 0
    while True:
        spline = sprot.readline()
        if not spline:
            break
        else:
            spline = spline.rstrip()
        spline = spline.split('\t')
        trline = trembl_dict.get(spline[1],'')
        if trline:
            trline[1][4] = spline[4]
            try:
                sp = int(spline[4])
            except:
                sp = 0
            try:
                tr = int(trline[1][5])
            except:
                tr = 0
            trline[1][6] = str(sp + tr)
        elif int(spline[4]) >= 1:
            spline.insert(5, '')
            spline.insert(6, spline[4])
            trembl_dict[spline[1]] = [spline[3].split()[0], spline]
        else:
            missing += 1
##            print(missing, 'Sprot taxon not in Trembl', spline[1], spline[4])
    fout = open(os.path.join(folder, 'uniprot_fasta_analyze.txt'), 'w')
    print(header, file = fout)
    line_list = list(trembl_dict.values())
    line_list.sort()
    for line in line_list:
        line = line[1]
        print('\t'.join(line), file=fout)
    [f.close() for f in [sprot, trembl, fout]]
    return

def get_taxon_from_name(db, name, dict_list):
    """Looks up ncbi taxonomy numbers by species name.
    Written by Phil Wilmarth, OHSU, 2009.
    """
    name2tax = dict_list[0]
    sci2tax = dict_list[1]
#    id2tax = dict_list[2]
#    name2id = dict_list[3]
    
    if db == 'nr':
        taxon = name2tax.get(name, -1)
        if taxon == -1:
            try:
                taxon = int(name.split('_')[2]) # this is for nr only
            except:
                pass
    else:
        taxon = sci2tax.get(name, name2tax.get(name, -1))
    return taxon

def time_stamp_logfile(message, file_obj):
    """Prints message and time stamp to a log file.
    Written by Phil Wilmarth, OHSU, 2009.
    """
    import time
    print('%s on %s' % (message, time.ctime()), file=file_obj)
    return

# end of module




