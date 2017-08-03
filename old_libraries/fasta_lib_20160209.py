"""'fasta_lib.py' Written by Phil Wilmarth, OHSU.
Copyright 2010, Oregon Health & Science University.
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
# 6/5/2009 PW - minor bug fixes for processing nr databases
#             - changed how reversed accessions are written for nr databases
# 7/8/2009 PW - revised fasta reader method to check amino acid characters
#             - added error checking flags to FasterReader & Protein methods
# 4/23/2010 PW - fixed findPeptide wrt I/L being indistinguishable
#              - improved adding decoy_string to accessions in reverseProtein
#              - added more support for Ref_Seq entries in NCBI databases
# 4/12/2011 PW - added a tryptic digest function
# 2011 PW - Switched to arrays to hold Gi to Taxon data to reduce memory use
# 2014 PW - Needed to increase size limit for reading arrays back in
# 10/2014 PW - Updated the FTP URLs for UniProt downloads
#
#
#================================
def taxon_cmd_line_checker(argv):
#================================
    """Checks command line arguments for correctness.
    Returns dictionary of taxon, name pairs or empty dictionary.

    Written by Phil Wilmarth, OHSU, 2009.
    """
    tax_dict = {}
    if argv[0].endswith('.py'):
        argv = argv[1:]
    #
    # need to have an even number of (integer, name) arguments
    #
    try:
        pairs = [(argv[i], argv[i+1]) for i in range(0, len(argv), 2)]
        for (taxon, name) in pairs:
            for char in taxon:
                if not char.isdigit():
                    break
            tax_dict[abs(int(taxon))] = name
    #
    # print usage information if error in format
    #
    except:
        print '\n   ### invalid command line argument format ###\n'
        print '   arguments must be a series of "taxonomy number" "name" pairs'
        print '   where "taxonomy number" is integer and "name" is text string.'
        print '   example: prompt>python program.py 9606 human 4932 yeast\n'
    #
    # return dictionary, it will be empty if any errors encountered
    #
    return tax_dict
#
#
#=================================
def string_cmd_line_checker(argv):
#=================================
    """Checks command line arguments for correctness.
    Returns dictionary of string, name pairs or empty dictionary.

    Written by Phil Wilmarth, OHSU, 2009.
    """
    str_dict = {}
    if argv[0].endswith('.py'):
        argv = argv[1:]
    #
    # need to have an even number of (string, name) arguments
    #
    try:
        pairs = [(argv[i], argv[i+1]) for i in range(0, len(argv), 2)]
        for (string, name) in pairs:
            str_dict[string] = name
    #
    # print usage information if error in format
    #
    except:
        print '\n   ### invalid command line argument format ###\n'
        print '   arguments must be a series of "string" "name" pairs'
        print '   where "string" is a text pattern and "name" is text string.'
        print '   Note: enclose "string" in double quotes if it has whitespace.'
        print '   example: prompt>python program.py "Homo sapiens" human\n'
    #
    # return dictionary, it will be empty if any errors encountered
    #
    return str_dict
#
#
#===================================================
def get_folder( default_location, title_string="" ):
#===================================================
    """Dialog box to browse to a folder.  Returns folder path.

    Usage: full_folder_name = get_folder(default_location, [title]),
        where "default_location" is a starting folder location,
        "title" is an optional message to list in the dialog box,
        and "full_folder_name" is the complete selected folder name.

    Written by Phil Wilmarth, 2008.
    """
    import os
    import Tkinter
    import tkFileDialog
    #
    # set up GUI elements
    #
    root = Tkinter.Tk()
    root.withdraw()
    try:
        root.tk.call('console', 'hide')
    except:
        pass
    #
    # set default title string if not passed
    #
    if title_string == "":   
        title_string = 'Select parent folder containing desired folders'
    #
    # create dialog box for folder selection
    #
    full_folder_name = tkFileDialog.askdirectory(parent=root, \
                                       initialdir=default_location, \
                                       title=title_string, mustexist=True)
    #
    # return full folder name
    #
    return full_folder_name   
    #
    # end
    #
#
#
#===============================================================
def get_file(default_location, extension_list, title_string=""):
#===============================================================
    """Dialog box to browse to a file.  Returns full file name.

    Usage: full_file_name = get_file(default_location, extension_list, [title]),
        where "default_location" is a starting folder location,
        extension_list is a list of (label, pattern) tuples,
        e.g. extension_list = [('Text files', '*.txt')],
        "title" is an optional message to list in the dialog box, and
        "full_file_name" is the complete name of the selected file.

    Written by Phil Wilmarth, OHSU, 2008.
    """
    import os
    import Tkinter
    import tkFileDialog
    #
    # set up GUI elements
    #
    root = Tkinter.Tk()
    root.withdraw()
    try:
        root.tk.call('console', 'hide')
    except:
        pass
    #
    # set default title string if not passed
    #
    if title_string == "":   
        title_string = 'Select a single FILE'
    #
    # create dialog box for file selection
    #
    filename = tkFileDialog.askopenfilename(parent=root, \
                                            initialdir=default_location, \
                                            filetypes=extension_list, \
                                            title=title_string)
    #
    # return full filename
    #
    return filename   
    #
    # end
    #
#
#
#================================================
def save_file(default_location, extension_list, \
              default_file='', title_string=""):
#================================================
    """Dialog box to save a file.  Returns full name of desired file.

    Usage: full_file_name = save_file(def_loc, ext_list, [def_file], [title]),
        where "def_loc" is a starting folder location,
        ext_list is a list of (label, pattern) tuples,
        e.g. extension_list = [('Text files', '*.txt')],
        "def_file" is an optional default filename,
        "title" is an optional message to list in dialog box, and
        "full_file_name" is the complete name of the desired file.

    Written by Phil Wilmarth, OHSU, 2009.
    """
    import os
    import Tkinter
    import tkFileDialog
    #
    # set up GUI elements
    #
    root = Tkinter.Tk()
    root.withdraw()
    try:
        root.tk.call('console', 'hide')
    except:
        pass
    #
    # set default title string if not passed
    #
    if title_string == "":   
        title_string = 'Select a single FILE'
    #
    # create dialog box for file selection
    #
    filename = tkFileDialog.asksaveasfilename(parent=root, \
                                              initialdir=default_location, \
                                              initialfile=default_file, \
                                              filetypes=extension_list, \
                                              title=title_string)
    #
    # return full filename
    #
    return filename   
    #
    # end
    #
#
#
#=============
class Protein:
#=============
    """Object to hold protein accession numbers, descriptions, and sequences.

    Methods:
        __init_:standard constructor, no parameters.
        readProtein: returns next protein from "fasta_reader"
        printProtein: prints sequence in FASTA format
        parseIPI: cleans up IPI entries
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
    #
    #------------------
    def __init__(self):
    #------------------
        """Basic constructor, no parameters.
        """
        #
        # bare bones __init__ function
        #
        self.accession = 'blank'
        self.new_acc = 'blank'
        self.description = 'blank'
        self.new_desc = 'blank'
        self.sequence = ''
        return
    #
    #-----------------------------------
    def readProtein(self, fasta_reader):
    #-----------------------------------
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
    #
    #------------------------------------------------
    def printProtein(self, file_obj=None, length=80):
    #------------------------------------------------
        """Prints FASTA protein entry to file (stdout is default).

        Usage: object.printProtein([file_obj=None, length=80]),
            where "object" is an instance of a Protein object, and
            "file_obj" is a file object (a value of None will print
            to standard out stream.  Optional "length" is number of
            characters per line for the protein sequence.

        Written by Phil Wilmarth, OHSU, 2009.
        """
        import sys
        if file_obj == None:
            file_obj = sys.stdout
        #
        # print new accession and new descriptor on first line
        #
        if self.new_desc == '':
            print >>file_obj, '>'+self.new_acc
        else:
            print >>file_obj, '>'+self.new_acc, self.new_desc
        #
        # initialize some things
        #
        char_count = 0
        char_line = ''
        #
        # build up sequence line with "length" characters per line
        #
        for char in self.sequence:
            if char_count < length:  # do not have "width" chars yet
                char_line += char
                char_count += 1
            else:                   # line is "width" long so print and reset
                print >>file_obj, char_line
                char_line = char
                char_count = 1
        #
        # print last sequence line (often less than "width" long) and return
        #
        if len(char_line):
            print >>file_obj, char_line
        #
        return
    #
    #--------------------------------------
    def parseIPI(self, KEEP_GENE_ID=False):
    #--------------------------------------
        """Parses IPI database accession numbers and descriptions.

        Written by Phil Wilmarth, OHSU, 2009.
        """
        #
        # extract the IPI number (keep version) and get rid of cross references
        #
        temp = self.accession.split('|')
        temp = [x.split(':')[1] for x in temp]
        self.new_acc = temp[0]
        other_acc = '|'.join(temp[1:])
        #
        # skip the taxonomy info since we already know species
        #
        try:
            new_desc = self.description.split(chr(01))[0]
            new_desc = new_desc.split('Tax_Id=')[1]
        except IndexError:
            pass
        temp = new_desc.split()
        new_desc = ' '.join(temp[1:])
        #
        # remove "Gene_Symbol" stuff if desired
        #
        new_desc = new_desc.replace('Gene_Symbol=', 'Gene=')
        if not KEEP_GENE_ID:    
            try:
                new_desc = new_desc.split('Gene=')[1]
            except IndexError:
                pass
            temp = new_desc.split()
            new_desc = ' '.join(temp[1:])
        #
        # add other accessions at end of description
        #
        if new_desc.endswith('.'):
            new_desc = new_desc[:-1]
        new_desc = '%s (%s).' % (new_desc, other_acc)
        self.new_desc = new_desc
        return
    #
    #---------------------------------------
    def parseNCBI(self, REF_SEQ_ONLY=False):
    #---------------------------------------
        """Parses NCBI nr database accession numbers and descriptions

        Written by Phil Wilmarth, OHSU, 2009.
        """
        #
        # keep the gi number (or RefSeq) for the accession, get rid of others
        #
        temp = self.accession.split('|')
        if REF_SEQ_ONLY:
            try:
                iacc = temp.index('ref')
            except ValueError:
                iacc = temp.index('gi')
        else:
            iacc = temp.index('gi')
        self.new_acc = '|'.join(temp[iacc:iacc+2])
        #
        # keep the first description string if more than one
        #
        new_desc = self.description.split(chr(01))[0]
        new_desc = new_desc.rstrip()
        #
        # add a period to the end of description
        #
        if not new_desc.endswith('.'):
            new_desc = new_desc + '.'
        self.new_desc = new_desc
        return
    #
    #---------------------------------------------
    def parseUniProt(self, KEEP_UNIPROT_ID=False):
    #---------------------------------------------
        """Parses UniProt database accession numbers and descriptions.

        Written by Phil Wilmarth, OHSU, 2009.
        """
        if len(self.accession.split('|')) < 3:
            print '   parseUniProt WARNING: fewer than 3 accession elements'
        #
        # if KEEP_UNIPROT_ID is set to True, keep the 
        # more human readible part of the accession numbers
        #
        if KEEP_UNIPROT_ID:
            self.new_acc = self.accession.split('|')[-1]
            acc_num = self.accession.split('|')[-2]
        else:
            self.new_acc = self.accession.split('|')[-2]
            acc_num = self.accession.split('|')[-1]
        #
        # keep the desciption text and get rid of the trailing stuff
        #
        new_desc = self.description.split(chr(01))[0]
        new_desc = new_desc.split('OS=')[0]
        new_desc = new_desc.rstrip()
        #
        # add other accession text to end of description
        #
        if new_desc.endswith('.'):
            new_desc = new_desc[:-1]
        self.new_desc = '%s (%s).' % (new_desc, acc_num)
        return
    #
    #-------------------
    def parseCONT(self):
    #-------------------
        """Parses contaminant accession numbers and descriptions

        Written by Phil Wilmarth, OHSU, 2009.
        """
        #
        # keep the contaminant tag for the accession, get rid of others
        #
        old_acc = ''
        temp = self.accession.split('|')
        if len(temp) > 1:
            self.new_acc = temp[0]
            old_acc = '|'.join(temp[1:])
        #
        # keep the first description string if more than one
        #
        new_desc = self.description.split(chr(01))[0]
        new_desc = new_desc.strip()
        #
        # add other accessions to description (if any)
        #
        if new_desc.endswith('.'):
            new_desc = new_desc[:-1]
        if old_acc:
            self.new_desc = '%s (%s).' % (new_desc, old_acc)
        else:
            self.new_desc = new_desc + '.'
        return
    #
    #--------------------------------------
    def reverseProtein(self, decoy_string):
    #--------------------------------------
        """Reverses protein sequence and returns new Protein object.

        Usage: rev_prot = object.reverseProtein(decoy_string),
            where "object" is a Protein object, "decoy_string" is the
            unique identifier text to add to the beginning of the 
            protein accesion number, and "rev_prot" is new Protein object.

        Written by Phil Wilmarth, OHSU, 2009.
        """
        import fasta_lib
        #
        # make sure decoy_string ends with an undescore
        #
        if not decoy_string.endswith('_'):
            decoy_string = decoy_string + '_'
        #
        # create a new Protein instance
        #
        rev_prot = Protein() 
        #
        # prefix the decoy_string to all important parts of accession
        #
        temp = self.accession.split('|')
        if self.accession.startswith('gi|'):    # modify nr DB accessions
            temp = [temp[i] for i in range(1, len(temp), 2)]
            temp = temp[0:1] # keep the gi number
        elif self.accession.startswith('sp|') or \
             self.accession.startswith('tr|'):  # modifiy Sprot or Trembl accs
            temp = temp[1:]
            temp = temp[0:1] # keep the accession number
        elif self.accession.startswith('IPI:'):     # modifiy IPI accessions
            temp = [x.split(':')[1] for x in temp]
            temp = temp[0:1] # keep the IPI number
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
        #
        # change the desciptions, too.
        #
        rev_prot.description = 'REVERSED.'
        rev_prot.new_desc = 'REVERSED.'
        #
        # now reverse the sequence
        #
        temp = list(self.sequence)
        temp.reverse()
        rev_prot.sequence = ''.join(temp)
        #
        return rev_prot
    #
    #--------------------------------------
    def molwtProtein(self, show_errs=True):
    #--------------------------------------
        """Returns protein molecular weight as the sum of average aa masses.
        If "show_errs" flag set, invalid amino acid characters are reported.

        Written by Phil Wilmarth, OHSU, 2009.
        """
        ave = { 'X':0.00000,
                'G':57.0519,
                'A':71.0788,
                'S':87.0782,
                'P':97.1167,
                'V':99.1326,
                'T':101.1051,
                'C':103.1388,
                'L':113.1594,
                'I':113.1594,
                'J':113.1594,
                'N':114.1038,
                'O':114.1472,
                'B':114.5962,
                'D':115.0886,
                'Q':128.1307,
                'K':128.1741,
                'Z':128.6231,
                'E':129.1155,
                'M':131.1926,
                'H':137.1411,
                'F':147.1766,
                'R':156.1875,
                'Y':163.1760,
                'W':186.2132,
                'U':168.053,
                '*':0.00000,
                '-':0.00000 }
        #
        # start with water and H+ masses, then add aa masses
        #
        bad_char = {}
        molwt = 18.01 + 1.007825
        for i in self.sequence:
            try:
                molwt += ave[i]
            except:     # keep track of bad characters
                bad_char[i] = True
        #
        bad_char = bad_char.keys()
        bad_char.sort()
        if len(bad_char) > 0 and show_errs:     # report bad chars if desired
            print '   WARNING: unknown symbol(s) (%s) in %s:\n%s' % \
                  (''.join(bad_char), self.accession, self.sequence)
        return molwt
    #
    #------------------------------------------
    def frequencyProtein(self, show_errs=True):
    #------------------------------------------
        """Returns aa frequency distrubution as a dictionary.
        If "show_errs" flag set, invalid amino acid characters are reported.
        
        Written by Phil Wilmarth, OHSU, 2009.
        """
        freq = {'X':0,
                'G':0,
                'A':0,
                'S':0,
                'P':0,
                'V':0,
                'T':0,
                'C':0,
                'L':0,
                'I':0,
                'J':0,
                'N':0,
                'O':0,
                'B':0,
                'D':0,
                'Q':0,
                'K':0,
                'Z':0,
                'E':0,
                'M':0,
                'H':0,
                'F':0,
                'R':0,
                'Y':0,
                'W':0,
                'U':0,
                '*':0,
                '-':0 }
        #
        # count the amino acids for all residues in sequence
        #
        bad_char = {}
        for i in self.sequence:
            try:
                freq[i] += 1
            except:     # keep track of bad characters
                bad_char[i] = True
        #
        bad_char = bad_char.keys()
        bad_char.sort()
        if len(bad_char) > 0 and show_errs: # report any bad chars, if desired
            print '   WARNING: unknown symbol(s) (%s) in %s:\n%s' % \
                  (''.join(bad_char), self.accession, self.sequence)
        return freq
    #
    #-----------------------
    def seqlenProtein(self):
    #-----------------------
        """Calculates protein sequence length.

        Written by Phil Wilmarth, OHSU, 2009.
        """
        return (len(self.sequence))
    #
    #------------------------------
    def findPeptide(self, peptide):
    #------------------------------
        """Calculates location of all 'peptide' matches in 'self.sequence.'

        Written by Phil Wilmarth, OHSU, 2009.
        """
        matches = []
        #
        # get rid of bounding residues, if any
        #
        try:
            peptide = peptide.split('.')[1]
        except IndexError:
            pass
        #
        # remove any modification symbols
        #
        base_pep = ''
        for c in peptide:
            if c.isalpha(): base_pep += c.upper()
        #
        # consider I/L indistinguishible
        #
        base_pep = base_pep.replace('I', 'j')
        base_pep = base_pep.replace('L', 'j')
        sequence = '-' + self.sequence + '-' # add bounding symbols
        sequence_fixed = sequence.replace('I', 'j')
        sequence_fixed = sequence_fixed.replace('L', 'j')
        #
        # find all matches of base_pep to protein sequence
        #
        if base_pep not in sequence_fixed:
            return []
        start = 0
        index = 0
        #
        while index != -1:
            index = sequence_fixed[start:].find(base_pep)
            start += index  # beginning residue number
            if index != -1: # skip this if no match
                end = start + len(base_pep) - 1 # ending residue no.
                #
                # add bounding aas, periods, and put back special chars
                #
                full_seq = sequence[start-1:start+len(base_pep)+1]
                full_seq = full_seq[0:1]+'.'+full_seq[1:-1]+'.'+full_seq[-1:]
                full_seq = full_seq.replace(base_pep, peptide)
                matches.append((start, end, full_seq))
                start = end + 1
        #
        # return the match list (empty list if no matches)
        #
        return matches
    #
    #------------------------------------
    def calcCoverage(self, peptide_list):
    #------------------------------------
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
    #
    # end class
    #
#
#
#=================
class FastaReader:
#=================
    """Reads FASTA entries from a file-like object.

    methods:
    __init__: basic constructor, no parameters.
    
    readProtein: reads one FASTA entry from a file object (text or zipped)
        arguments are "next_protein" and "file_obj"
        returns True (next protein) or False (EOF or not FASTA).

    written by Phil Wilmarth, OHSU, 2009.
    """
    #
    #------------------------------
    def __init__(self, fasta_file):
    #------------------------------
        """Basic constructor function.  No parameters

        self._last_line retains the previous '>' line and
        self._valid is a dictionary of valid protein FASTA chars.
        """
        #
        import os
        import gzip
        import fasta_lib
        #
        # attribute to save last line from previous read
        #
        self._last_line = 'start value'
        self._file_obj = None
        self._fasta_file = fasta_file
        #
        # list of valid amino acid characters
        #
        self._valid = {'X':True, 'G':True, 'A':True, 'S':True, 'P':True,\
                      'V':True, 'T':True, 'C':True, 'L':True, 'I':True,\
                      'J':True, 'N':True, 'O':True, 'B':True, 'D':True,\
                      'Q':True, 'K':True, 'Z':True, 'E':True, 'M':True,\
                      'H':True, 'F':True, 'R':True, 'Y':True, 'W':True,\
                      'U':True, '*':True, '-':True }
        #
        # get file object and save as attribute
        #
        if not os.path.exists(fasta_file):
            fasta_file = fasta_lib.get_file(default_location, \
                                              extension_list, \
                                              title_string="Select FASTA file")
        try:
            if fasta_file.endswith('.gz'):
                self._file_obj = gzip.open(fasta_file, 'rb')
            else :
                self._file_obj = open(fasta_file, 'Ur')
        except IOError:
            print '   WARNING: Fasta database could not be opened!'
            raise
        return
    #
    #-------------------------------------------------------------
    def readNextProtein(self, next_protein, check_for_errs=False):
    #-------------------------------------------------------------
        """Loads one FASTA protein text entry into a Protein object.

        Returns True (protein entry found) or False (end of file).
        If "check_for_errs" flag is set, amino acid chars are checked.

        Written by Phil Wilmarth, OHSU, 2009.
        """
        #
        # at first call, start reading lines
        #
        if self._last_line == 'start value':
            self._last_line = self._file_obj.readline()
            if not self._last_line:
                self._file_obj.close()
                return(False)
            self._last_line = self._last_line.strip()
        #
        # get next protein's info from _last_line
        #
        if self._last_line.startswith('>'):
            next_protein.accession = self._last_line.split()[0][1:]
            next_protein.new_acc = next_protein.accession
            start = len(next_protein.accession)+2
            next_protein.description = self._last_line[start:]
            next_protein.new_desc = next_protein.description
        #
        # return if empty line (EOF) or non-description line
        #
        else:
            self._file_obj.close()
            return(False)                    
        #
        # reset variables and read in next entry
        #
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
            #
            # stop reading at next descriptor line (and save line)
            #
            if line.startswith('>'):
                self._last_line = line.strip()
                #
                # report bad characters if conditions were met
                #
                bad_char = bad_char.keys()
                bad_char.sort()
                if len(bad_char) > 0 and check_for_errs:
                    print '   WARNING: unknown symbol(s) (%s) in %s' % \
                          (''.join(bad_char), next_protein.accession)
                break
            #
            # add next sequence line to protein's sequence
            #
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
        #
        # return (protein info retained in next_protein)
        #
        return(True)
        #
        # end
        #
    #
    # end class
    #
#
#
#=========================
def get_uniprot_version():
#=========================
    """Gets UniProt version numbers from online release notes.

    Written by Phil Wilmarth, OHSU, 2009.
    FTP addresses updated -PW 2014.
    """
    import urllib
    #
    # set up to read the online UniProt release notes file
    #
    print '...getting database version numbers...'
    versions = {'uniprot':'XX.X', 'sprot':'XX.X', 'trembl':'XX.X'}
    address = 'ftp://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/complete/reldate.txt'
    reldate = urllib.urlopen(address)
    #
    # read release notes file and get UniProt, etc. version numbers
    #
    for line in reldate.readlines():
        if 'UniProt Knowledgebase' in line:
            versions['uniprot'] = line.split()[3].replace('_', '.')
        elif 'Swiss-Prot' in line:
            versions['sprot'] = line.split()[2].replace('_', '.')
        elif 'TrEMBL' in line:
            versions['trembl'] = line.split()[2].replace('_', '.')
    #
    return(versions)
#
#
#======================
def get_ipi_versions():
#======================
    """Gets IPI version numbers from online release notes.

    Written by Phil Wilmarth, OHSU, 2009.
    """
    import urllib
    #
    # set up to read the online UniProt release notes file
    #
    print '...getting IPI database version numbers...'
    versions = {}
    entries = {}
    address = 'ftp://ftp.ebi.ac.uk/pub/databases/IPI/current/README'
    readme = urllib.urlopen(address)
    #
    # read release readme file and get version numbers
    #
    species_dict = {'human':'HUMAN', 'mouse':'MOUSE', 'rat':'RAT',
                    'zebrafish':'DANRE', 'arabidopsis':'ARATH',
                    'chicken':'CHICK', 'cow':'BOVIN'}
    parse = False
    for line in readme.readlines():
        if line.startswith('Current release'):
            parse = True
        if parse:
            line = line.rstrip()
            temp = line.split()
            if len(temp) >= 3:
                if temp[-3] in species_dict.keys():
                    versions[species_dict[temp[-3]]] = temp[-2]
                    entries[species_dict[temp[-3]]] = temp[-1]
        if line.startswith('This release was built using'):
            break
    #
    return(versions, entries)
#
#
#======================================
def download_ipi(db, folder, versions):
#======================================
    """Downloads (if necessary) IPI fasta file to "folder".

    Written by Phil Wilmarth, OHSU, 2009.
    """
    import os
    import urllib
    import socket
    #
    # check if files are already downloaded, if not fetch them from ipi site
    #
    socket.setdefaulttimeout(120.)
    print '...downloading IPI databases and release notes...'
    db_name = 'ipi.%s.%s.fasta.gz' % (db.upper(), versions[db])
    base_address = 'ftp://ftp.ebi.ac.uk/pub/databases/IPI/current/'
    db_address = 'ipi.%s.fasta.gz' % (db.upper(),)
    db_address = base_address + db_address
    files_addresses = [ (os.path.join(folder, db_name), db_address),
                         (os.path.join(folder, 'README'),
                         'ftp://ftp.ebi.ac.uk/pub/databases/IPI/current/README')]
    files_addresses.reverse()
    for (file_name, address) in files_addresses:
        if os.path.exists(file_name):
            print '...existing version is current'
            pass
        else:
            print '...downloading', file_name
            x = reporter()
            try:
                urllib.urlretrieve(address, file_name, reporthook=x.report)
            except:
                print '...WARNING: download may have hung at EOF or other error'
            finally:
                urllib.urlcleanup()
    return
#
#
#==========================================
def download_uniprot(db, folder, versions):
#==========================================
    """Downloads (if necessary) UniProt FASTA, taxon files to "folder".

    Written by Phil Wilmarth, OHSU, 2009.
    FTP addresses updated -PW 2014.
    """
    import os
    import urllib
    import socket
    #
    # check if files are already downloaded, if not fetch them from uniprot site
    #
    socket.setdefaulttimeout(120.)
    print '...downloading databases and taxonomy files...'
    db_name = 'uniprot_%s_%s.fasta.gz' % (db, versions[db],)
    base_address = 'ftp://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/complete/'
    db_address = 'uniprot_%s.fasta.gz' % (db,)
    db_address = base_address + db_address
    files_addresses = [ (os.path.join(folder, db_name), db_address),
                         (os.path.join(folder, 'taxdump.tar.gz'),
                         'ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz'),
                         (os.path.join(folder, 'speclist.txt'),
                         'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/complete/docs/speclist.txt'),
                         (os.path.join(folder, 'relnotes.htm'),
                        'ftp://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/complete/docs/relnotes.htm') ]
    files_addresses.reverse()
    for (file_name, address) in files_addresses:
        if os.path.exists(file_name):
            pass
        else:
            print '...downloading', file_name
            x = reporter()
            try:
                urllib.urlretrieve(address, file_name, reporthook=x.report)
            except:
                print '...WARNING: download may have hung at EOF or other error'
            finally:
                urllib.urlcleanup()
    return
#
#
#============================
def download_ncbi(nr_folder):
#============================
    """Downloads (if necessary) the ncbi FASTA and taxon files to "nr_folder".

    Written by Phil Wilmarth, OHSU, 2009.
    """
    import os
    import urllib
    import socket
    #
    # check if files are already downloaded, if not fetch them from ncbi site
    #
    nr_name = os.path.split(nr_folder)[1] + '.gz'
    if os.path.exists(os.path.join(nr_folder, 'nr.gz')):
        os.rename(os.path.join(nr_folder, 'nr.gz'), \
                  os.path.join(nr_folder, nr_name))
    files_addresses = [ (os.path.join(nr_folder, nr_name),
                    'ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nr.gz'),
                   (os.path.join(nr_folder, 'gi_taxid_prot.dmp.gz'),
                    'ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_prot.dmp.gz'),
                   (os.path.join(nr_folder, 'taxdump.tar.gz'),
                    'ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz') ]
    files_addresses.reverse()
    socket.setdefaulttimeout(120.)
    for (file_name, address) in files_addresses:
        if os.path.exists(file_name):
            pass
        else:
            print '...downloading', file_name
            x = reporter()
            try:
                urllib.urlretrieve(address, file_name, reporthook=x.report)
            except:
                print '...WARNING: download may have hung at EOF or other error'
            finally:
                urllib.urlcleanup()
    return
#
#----------------
class reporter():
#----------------
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
            print '......%s of %s bytes (%.2f%%)' % \
                  (sub_total, size, float(100*sub_total)/size)
        return
    #
    # end class
    #
#
#
#============================================================
def expand_species(folder, db, taxon_dict, \
                   min_sequence_count, min_seq_per_species, \
                   REF_SEQ_ONLY=False):
#============================================================
    """Expands any taxon nodes numbers into all member taxon numbers.

    Written by Phil Wilmarth, OHSU, 2009.
    """
    import os
    import tarfile
    #
    VERBOSE = False
    #
    # open taxonomy nodes file
    #
    print '...making taxonomy nodes dictionary...'
    archive_name = os.path.join(folder, 'taxdump.tar.gz')
    archive = tarfile.open(archive_name)
    nodes = archive.extractfile('nodes.dmp')
    #
    # read nodes file and save child to parent taxon mappings
    #
    child_to_parent = {}
    while True:
        line = nodes.readline()
        if not line:
            break
        else:
            line = line.rstrip()
        item = line.split('\t|\t')
        child_to_parent[int(item[0])] = int(item[1])    
    nodes.close()
    #
    # open the fasta_analysis.txt file
    #
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
    #
    # save species taxons and sequence counts
    #
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
    #
    # see if we have any group taxon numbers
    #
    group_expand = {}
    for alltax in [x for x in species_counts.keys() \
                   if species_counts[x] >= min_seq_per_species]:
        tree = []
        tree.append(alltax)
        parent = alltax
        while parent != 1:
            try:
                parent = child_to_parent[parent]
            except KeyError:
                if VERBOSE:
                    print '...WARNING: lookup of taxon %s failed' % \
                          (parent,)
                break
            tree.append(parent)
        #
        # see if any of our taxon dictionary numbers are in lineage
        #
        for tax in taxon_dict.keys():
            if (tax in tree) and (alltax != tax):
                    taxon_dict[alltax] = taxon_dict[tax]
                    try:
                        group_expand[tax].append(alltax)
                    except KeyError:
                        group_expand[tax] = [alltax]
    #
    # report any expanded taxononmy numbers
    #
    for tax in taxon_dict.keys():
        if len(group_expand.get(tax, [])) >= 1:
            print '...taxon %s (%s direct seqs) was a node with children:' \
                  % (tax, species_counts.get(tax, 0))
            for child in group_expand[tax]:
                print '......taxon %s has %s sequences' % \
                      (child, species_counts[child])                
    #
    # check taxon dictionary and remove taxons with too few sequences
    #
    for tax in taxon_dict.keys():
        if species_counts.get(tax, 0) <= min_sequence_count:
            print '...WARNING: taxon number %s had too few proteins' % (tax,)
            del taxon_dict[tax]
    return
#
#
#=============================================
def tax_lookup_by_gi(gi, gi_array, tax_array):
#=============================================
    """Binary lookup of taxon number given a gi number.

    Written by Phil Wilmarth, OHSU, 2009.
    Adapted algorithm from John Zelle "Python Programming:
    An Introduction to Computer Science", Ch. 13.1.3,
    Franklin, Beedle & Associates, Wilsonville OR, 2004.
    """
    low = 0
    high = len(gi_array) - 1
    while low <= high:
        mid = (low + high) / 2
        if gi_array[mid] == gi:
            return tax_array[mid]
        elif gi < gi_array[mid]:
            high = mid - 1
        else:
            low = mid + 1
    return -1
#
#
#=================
class GiToTaxon():
#=================
    """Object to map NCBI gi numbers to taxonomy numbers.

    Methods:
        __init__: loads arrays from taxonomy files (or reloads)
        get(gi, default): return taxon number of "gi" or "default"

    Written by Phil Wilmarth, OHSU, 2009.
    """
    #
    #-----------------------------
    def __init__(self, nr_folder):
    #-----------------------------
        import array
        import os
        import gzip
        #
        # arrays to hold gi numbers and taxon numbers
        #
        self.gi_array = array.array('i')
        self.tax_array = array.array('i')
        return
    #
    #-----------------------------------
    def create_or_load(self, nr_folder):
    #-----------------------------------
        """Creates or reloads a GiToTaxon object.
        """
        import array
        import os
        import gzip
        import cPickle
        #
        # build filenames
        #
        nr_name = os.path.split(nr_folder)[1] + '.gz'
        nr_file = os.path.join(nr_folder, nr_name)
        fasta_file = gzip.open(nr_file,'rb')
        gi_file = os.path.join(nr_folder, 'gi_array_binary.txt')
        tax_file = os.path.join(nr_folder, 'tax_array_binary.txt')
        pickle_file = os.path.join(nr_folder, 'gi_to_taxon_pickle.txt')
        #
        # check for array file existance
        #
        if os.path.exists(gi_file) and os.path.exists(tax_file):
            #
            # if arrays have been saved, then reload from disk
            #
            print '...reloading gi_to_taxon mapping arrays...'
            try:
                f = open(gi_file, 'rb')
                self.gi_array.fromfile(f, 999000000)
            except EOFError:
                f.close()
            try:
                f = open(tax_file, 'rb')
                self.tax_array.fromfile(f, 999000000)
            except EOFError:
                f.close()
            print '...%s gi and taxon numbers read from disk' % len(self.gi_array)
        elif os.path.exists(pickle_file):
            print '...reloading gi_to_taxon pickled file...'
            f = open(pickle_file, 'rb')
            p = cPickle.Unpickler(f)
            temp = p.load()
            f.close()
            self.gi_array = temp.gi_array
            self.tax_array = temp.tax_array
        else:
            print '...making gi_to_taxon mapping information...'
            #
            # read the gi number to taxonomy number data into arrays
            #
            fname = os.path.join(nr_folder, 'gi_taxid_prot.dmp.gz')
            fin = gzip.open(fname,'rb')
            line_count = 0
            while True:
                line = fin.readline()
                if not line:
                    break
                else:
                    line = line.rstrip()
                line_count += 1
                if (line_count % 1000000) == 0:
                    print '......(%s gi_to_taxon lines read)' % (line_count,)
                gi = int(line.split('\t')[0])
                tax = int(line.split('\t')[1])
                self.gi_array.append(gi)
                self.tax_array.append(tax)
            fin.close()
            f = open(gi_file, 'wb')
            self.gi_array.tofile(f)
            f.close()
            f = open(tax_file, 'wb')
            self.tax_array.tofile(f)
            f.close()
            print '...%s gi and taxon number saved to disk' % len(self.gi_array)
        #
        # make sure arrays are same size
        #
#        print 'lenght of arrrays', len(self.gi_array)
        if len(self.gi_array) != len(self.tax_array):
            print '   WARNING: gi array and tax array are not in sync!'
        return
    #
    #--------------------------
    def get(self, gi, default):
    #--------------------------
        """Binary lookup of taxon number given a gi number.

        Written by Phil Wilmarth, OHSU, 2009.
        Adapted algorithm from John Zelle "Python Programming:
        An Introduction to Computer Science", Ch. 13.1.3,
        Franklin, Beedle & Associates, Wilsonville OR, 2004.
        """
        low = 0
        high = len(self.gi_array) - 1
        while low <= high:
            mid = (low + high) / 2
            if self.gi_array[mid] == gi:
                return self.tax_array[mid]
            elif gi < self.gi_array[mid]:
                high = mid - 1
            else:
                low = mid + 1
        #
        return default
    #
    # end class
    #
###
###
###===============================
##def make_gi_to_taxon(nr_folder):
###===============================
##    """Makes and returns the gi_to_taxon dictionary.
##    Also pickles the dictionary for later use (slow to read nr).
##
##    Written by Phil Wilmarth, OHSU, 2009.
##    """
##    import os
##    import gzip
##    import cPickle
##    import array
##    import fasta_lib
##    #
##    # load dictionary if it exists to save time
##    #
##    gi_to_taxon = GiToTaxon()
##    gi_to_taxon.create_or_load(nr_folder)
##    
##    save_name = os.path.join(nr_folder, 'gi_to_taxon_dict_pickle.txt')
##    if os.path.exists(save_name):
##        print '...reading gi_to_taxon dictionary from disk...'
##        f = open(save_name, 'rb')
##        p = cPickle.Unpickler(f)
##        gi_to_taxon = p.load()
##        return (gi_to_taxon, 0)
##    #
##    # else create the dictionary from scratch
##    #
##    print '...making gi_to_taxon dictionary...'
##    nr_name = os.path.split(nr_folder)[1] + '.gz'
##    nr_file = os.path.join(nr_folder, nr_name)
##    fasta_file = gzip.open(nr_file,'rb')
##    #
##    # load the gi number to taxonomy number data into arrays
##    #
##    gi_array = array.array('i')
##    tax_array = array.array('i')
##    fname = os.path.join(nr_folder, 'gi_taxid_prot.dmp.gz')
##    fin = gzip.open(fname, 'rb')
##    while True:
##        line = fin.readline()
##        if not line:
##            break
##        else:
##            line = line.rstrip()
##        gi = int(line.split('\t')[0])
##        tax = int(line.split('\t')[1])
##        gi_array.append(gi)
##        tax_array.append(tax)
##    fin.close()
##    if len(gi_array) != len(tax_array):
##        print '   WARNING: gi array and tax array are not in sync!'
##    #
##    # start reading nr.fasta (zipped) file, look for accessions
##    #
##    print '...reading nr database'
##    prot = 0
##    gi_to_taxon = {}
##    while True:
##        line = fasta_file.readline()
##        if not line:
##            break
##        else:
##            line = line.rstrip()
##        if line.startswith('>'):
##            prot += 1
##            if (prot % 100000) == 0:
##                print '...(%s proteins read so far)' % (prot,)
##            tax_list = []
##            for header in line.split(chr(01)):
##                gi = header.split()[0]
##                gi = int(gi.split('|')[1])
##                tax = fasta_lib.tax_lookup_by_gi(gi, gi_array, tax_array)
##                if tax in tax_list:
##                    continue
##                else:
##                    tax_list.append(tax)
##                    gi_to_taxon[gi] = tax        
##    fasta_file.close()
##    #
##    # save the dictionary to disk
##    #
##    save_name = os.path.join(nr_folder, 'gi_to_taxon_dict_pickle.txt')
##    save_file = open(save_name, 'wb')
##    cPickle.dump(gi_to_taxon, save_file, protocol=2)
##    save_file.close()
##    #
##    # return dictionary, nr protein counter
##    #
##    return (gi_to_taxon, prot)
#
#
#==================================
def make_taxon_to_sci_name(folder):
#==================================
    """Makes the taxon_to_sci_name dictionary.

    Written by Phil Wilmarth, OHSU, 2009.
    """
    import tarfile
    import os
    #
    print '...making taxon_to_sci_name dictionary...'
    archive_name = os.path.join(folder, 'taxdump.tar.gz')
    archive = tarfile.open(archive_name)
    names = archive.extractfile('names.dmp')
    taxon_to_name = {}
    #
    # read file and save names from 'scientific name' lines
    #
    line = names.readline().rstrip()
    while line:
        item = line.split('\t')
        if item[6] == 'scientific name':
            taxon_to_name[int(item[0])] = item[2]
        line = names.readline().rstrip()
    names.close()
    #
    # there may be some gi numbers that have taxon id of zero
    #
    taxon_to_name[0] = 'Zero_taxon_number'
    #
    return taxon_to_name
#
#
#=================================
def make_uniprot_to_taxon(folder):
#=================================
    """Makes sci_to_taxon and id_to_taxon dictionaries from "speclist.txt".

    Written by Phil Wilmarth, OHSU, 2009.
    """
    import os
    #
    print '...making scientific names and IDs to taxon dictionaries...'
    speclist = open(os.path.join(folder, 'speclist.txt'), 'r')
    sci_to_taxon = {}
    id_to_taxon = {}
    #
    # read file and save names from 'scientific name' lines
    #
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
    #
    # close file and return dictionaries
    #
    speclist.close()
    return(sci_to_taxon, id_to_taxon)
#
#
#===================================
def make_all_names_to_taxon(folder):
#===================================
    """Makes the all_names_to_taxon dictionary.

    Written by Phil Wilmarth, OHSU, 2009.
    """
    import tarfile
    import os
    #
    print '...making all_names_to_taxon dictionary...'
    archive_name = os.path.join(folder, 'taxdump.tar.gz')
    archive = tarfile.open(archive_name)
    names = archive.extractfile('names.dmp')
    all_names_to_taxon = {}
    #
    # read file and save taxonomy numbers for all names
    #
    while True:
        line = names.readline().rstrip()
        if not line: break          
        item = line.split('\t')
        name = item[2].replace('"','')
        name = name.lstrip()
        all_names_to_taxon[name] = int(item[0])
    #
    # close archive and return dictionary
    #
    names.close()
    return(all_names_to_taxon)
#
#
#============================================
def uniprot_species_frequency(database_name):
#============================================
    """Compiles species frequency info from Sprot or Trembl databases.

    Written by Phil Wilmarth, OHSU, 2009.
    """
    import gzip
    import fasta_lib
    #
    # read all of the protein descriptions and parse out species names
    #
    if 'sprot' in database_name:
        print '...scanning sprot database for species names...'
    else:
        print '...scanning trembl database for species names (slow)...'
    name_freq = {}
    name_to_spec_id = {}
    prot_count = 0
    f = gzip.open(database_name, 'rb')
    while True:
        line = f.readline()
        if not line: break
        #
        # get species name, id; save in dictionary; make frequency totals
        #
        if line.startswith('>'):
            prot_count += 1
            if (prot_count % 500000) == 0:
                print '......(%s proteins read...)' % (prot_count,)
            (spec_id, name) = uniprot_parse_line(line)
            name_to_spec_id[name] = spec_id
            fasta_lib.add_or_increment(name, name_freq)            
    #
    return(name_freq, name_to_spec_id, prot_count)
#
#
#============================
def uniprot_parse_line(line):
#============================
    """Parses UniProt description lines for species IDs and names.

    Written by Phil Wilmarth, OHSU, 2009.
    """
    name = 'Unknown species'
    spec_id = '_Z1234'
    #
    # get the species ID
    #
    spec_id = line.split()[0]
    spec_id = spec_id.split('|')[-1]
    spec_id = spec_id.split('_')[1]
    spec_id = '_' + spec_id
    #
    # get the species name if present
    #
    if len(line.split('OS=')) == 1:
        pass
    else:
        name = line.split('OS=')[1]
        if len(name.split('=')) == 1:
            name = name.rstrip()
        else:
            name = name.split('=')[0][:-2]
            name = name.rstrip()
    #
    return(spec_id, name)            
#
#
#================================================
def add_or_increment(species_name, species_dict):
#================================================
    """Increments species count (if existing) or adds species to dictionary.

    Written by Phil Wilmarth, OHSU, 2009.
    """
    #
    # add new species name or increment counter if name already exists
    #
    try:
        species_dict[species_name] += 1
    except:
        species_dict[species_name] = 1
    return
#
#
#==============================
def sort_species(species_dict):
#==============================
    """Sorts a dictionary of species frequencies by genus name then frequency.

    Written by Phil Wilmarth, OHSU, 2009.
    """
    import operator
    #
    # first sort by Genus name using the "sort" attribute
    #
    sorted_list = species_dict.items()
    sorted_list.sort()
    return sorted_list
##    #
##    # sort the Specific names by inverse protein sequence count
##    # we need to define the Genus subsections of the database
##    #
##    groups = []
##    current = sorted_list[0][0]
##    start = 0
##    end = 0
##    #
##    # find out where each Genus starts and stops in the list
##    #
##    for i, species in enumerate(sorted_list):
##        end += 1
##        if species[0] != current:
##            groups.append((start, end-1))
##            start = end-1
##            current = species[0]
##    groups.append((start, end))
##    #
##    # now we can sort by sequence count ("second_sort" attribute)
##    # and then reverse for descending order
##    #
##    for (start, end) in groups:
##        temp = [(y, x) for (x, y) in sorted_list[start:end]]
##        temp.sort()
##        temp.reverse()
##        sorted_list[start:end] = [(y, x) for (x, y) in temp]
##    #
##    return sorted_list
#
#
#====================
def virus_test(name):
#====================
    """Tests whether "virus" is in a species name or not.

    Written by Phil Wilmarth, OHSU, 2009.
    """
    if 'virus' in name:
        return 'TRUE'
    else:
        return 'FALSE'
    # end
#
#
#======================================================
def save_species_info_nr(folder, name_freq, name2tax, \
                         ref2freq=None, minimum=0):
#======================================================
    """Writes taxonomy and species name info to file.

    Written by Phil Wilmarth, OHSU, 2009.
    """
    import os
    import fasta_lib
    #
    # sort the species names and write to file
    #
    print '...writing species analysis results file...'
    sorted_list = fasta_lib.sort_species(name_freq)
    fout_name = os.path.join(folder, 'nr_fasta_analyze.txt')
    fout = open(fout_name, 'w')
    sformat = '%s\t%s\t%s\t%s\t%s\t%s\t%s'
    print >>fout, sformat % ('=SUBTOTAL(109,A2:A65000)', 'Taxon_ID', \
                             'Species_Name', 'Sequence_Count', \
                             'RefSeq_Count', 'Word_Count', 'Is_virus')
    dict_list = [name2tax, None, None, None]
    for name, count in sorted_list:
        if int(count) >= minimum:
            taxon = fasta_lib.get_taxon_from_name('nr', name, dict_list)
            #
            ref_count = ref2freq.get(taxon, '')
            print >>fout, sformat % (str(1), taxon, name, count,\
                                     ref_count, len(name.split()), \
                                     fasta_lib.virus_test(name))
    fout.close()
    return
#
#=====================================================================
def save_species_info(db, folder, name_freq, name2tax, sci2tax=None, \
                      id2tax=None, name2id=None, minimum=0):
#=====================================================================
    """Writes taxonomy and species name info to file.

    Written by Phil Wilmarth, OHSU, 2009.
    """
    import os
    import fasta_lib
    #
    # sort the species names and write to file
    #
    print '...writing species analysis results file...'
    sorted_list = fasta_lib.sort_species(name_freq)
    fout_name = os.path.join(folder, db+'_fasta_analyze.txt')
    fout = open(fout_name, 'w')
    sformat = '%s\t%s\t%s\t%s\t%s\t%s\t%s'
    print >>fout, sformat % ('=SUBTOTAL(109,A2:A65000)', 'Taxon_ID', \
                             'Species_ID', 'Species_Name', 'Sequence_Count', \
                             'Word_Count', 'Is_virus')
    dict_list = [name2tax, sci2tax, id2tax, name2id]
    for name, count in sorted_list:
        if int(count) >= minimum:
            taxon = fasta_lib.get_taxon_from_name(db, name, dict_list)
            #
            try:
                ID = name2id[name]
            except:
                ID = ' '
            print >>fout, sformat % (str(1), taxon, ID, name, count,\
                                     len(name.split()), \
                                     fasta_lib.virus_test(name))                
    fout.close()
    return
#
#
#==================================
def combine_analysis_files(folder):
#==================================
    """Program to add Sprot columns to Trembl analysis file

    Written by Phil Wilmarth, OHSU, 2009.
    """
    #
    import os
    trembl = open(os.path.join(folder, 'trembl_fasta_analyze.txt'), 'r')
    sprot = open(os.path.join(folder, 'sprot_fasta_analyze.txt'), 'r')
    header = trembl.readline().rstrip()
    new_header = header.split('\t')
    new_header[4] = 'Trembl_Count'
    new_header.insert(4, 'Sprot_Count')
    new_header.insert(6, 'Total_Count')
    header = '\t'.join(new_header)
    #
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
    #
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
    #        print missing, 'Sprot taxon not in Trembl', spline[1], spline[4]
    #
    out = open(os.path.join(folder, 'uniprot_fasta_analyze.txt'), 'w')
    print >>out, header
    line_list = trembl_dict.values()
    line_list.sort()
    for line in line_list:
        line = line[1]
        print >>out, '\t'.join(line)
    return
#
#
#============================================
def get_taxon_from_name(db, name, dict_list):
#============================================
    """Looks up ncbi taxonomy numbers by species name.

    Written by Phil Wilmarth, OHSU, 2009.
    """
    name2tax = dict_list[0]
    sci2tax = dict_list[1]
    id2tax = dict_list[2]
    name2id = dict_list[3]
    #
    if db == 'nr':
        taxon = name2tax.get(name, -1)
        if taxon == -1:
            try:
                taxon = int(name.split('_')[2]) # this is for nr only
            except:
                pass
    else:
        taxon = sci2tax.get(name, name2tax.get(name, -1))
    return(taxon)
#
#
#=========================================
def time_stamp_logfile(message, file_obj):
#=========================================
    """Prints message and time stamp to a log file.

    Written by Phil Wilmarth, OHSU, 2009.
    """
    import time
    print >>file_obj, '%s on %s' % (message, time.ctime())
    return
#
# end of module
#




