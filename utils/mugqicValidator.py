#!/usr/bin/env python

################################################################################
# Copyright (C) 2014, 2022 GenAP, McGill University and Genome Quebec Innovation Centre
#
# This file is part of MUGQIC Pipelines.
#
# MUGQIC Pipelines is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MUGQIC Pipelines is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with MUGQIC Pipelines.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

import getopt
import os
import sys
import re


def help():
        print('''
        ---------------------------------------------------------------------------
                                            HELP
        ---------------------------------------------------------------------------
        
        mugqicValidator.py validates the basic structure and integrity of files used
        by the GenAP pipelines. User must provide a readset file or a design file\n
        -r     --readset     myReadsetFile
        -d     --design      myDesignFile
        -h     --help



        usage: mugqicValidator.py -r myReadsetFile
               mugqicValidator.py -d myDesignFile

        ''')


def error():
        print('''
        Error---------Error---------Error---------Error---------Error---------Error
        ''')


def loadFile(myFile):
    '''
    loads file
    '''
    s = []
    f = open(myFile)
    for line in f:
        s.append(line.split("\t"))
    return s


def unicodeCheck(stringArray):
    '''
    Checks if there are non-ascii characters in file
    '''
    Pass = True
    for row in stringArray:
        try:
            [s.encode(encoding='utf-8').decode('ascii') for s in row]
        except UnicodeDecodeError:
            Pass = False
            error()
            print("File contains non English characters in line:\n" + str(row) + "at element:\n" + str(s))
    return Pass


def rowSizeCheck(stringArray):
    '''
    Checks that all rows in the file hasve the same size
    '''
    columns = []
    for row in stringArray:
        columns.append(len(row))
    if (all([columns[0] == n for n in columns])):
        return True
    else:
        error()
        print("Rows are not all the same size! please make all the rows of equal size!")
        return False
 
def trailingSpacesCheck(stringArray):
    '''
    Checks that words don't have spaces around them
    '''
    Pass = True
    for row in stringArray:
        for s in row:
            if re.search(r"[ \t]", s):
                error()
                print("There are trailing spaces around the following word: " + s + " ; Please fix it and try again!")
                Pass = False
    return Pass


def DosNewline(stringArray):
    '''
    Checks if the file has carriage returns and Dos newline characters
    '''
    Pass = True
    for row in stringArray:
        for s in row:
            if re.search(r"[\r]", s):
                error()
                print("There is a carriage return after the following word: " + s + " ; To fix it, you can run 'dos2unix myfile'. Please fix it and try again!")
                Pass = False
    return Pass


def readsetHeaderCheck(h):
    '''
    Checks that the readset file header is correct
    '''
    header = [
        ['Sample', 'Readset', 'Library', 'RunType', 'Run', 'Lane', 'Adapter1', 'Adapter2', 'QualityOffset', 'BED', 'FASTQ1', 'FASTQ2', 'BAM\n'],
        ['Sample', 'Readset', 'MarkName', 'MarkType', 'Library', 'RunType', 'Run', 'Lane', 'Adapter1', 'Adapter2', 'QualityOffset', 'BED', 'FASTQ1', 'FASTQ2', 'BAM\n']
    ]
    if h in header:
        return True
    else:
        error()
        print("File header is not correct. Please revise and try again!")
        return False

        
def readsetRunTypeCheck(stringArray):
    '''
    Checks if Run type is "SINGLE_END" or "PAIRED_END"
    '''
    RunType = [row[3] for row in stringArray][1:]
    if (all([(runtype == "SINGLE_END" or runtype == "PAIRED_END") for runtype in RunType])):
        return True
    else:
        error()
        print("RunType is not equal to 'SINGLE_END' or 'PAIRED_END' for all rows. Please correct before running pipelines!")
        return False


def readsetAdapterCheck(stringArray):
    '''
    Checks if Adapter sequences are a combination of {A,C,T,G}
    '''
    Pass = True
    entryOps = ['A','C','T','G']
    for row in stringArray[1:]:
        adapters = [s.strip().upper() for s in row[6:8]]
        for adapter in adapters:
            if not all([nucleotide in entryOps for nucleotide in adapter]):
                Pass = False
                error()
                print(str(adapter) + " contains values that are not permitted (" + str(entryOps) + "). Please fix this and try again!")
    return Pass



def readsetIntegerCheck(stringArray):
    '''
    Checks if QualityOffset is an integer
    '''
    Pass = True
    for row in stringArray[1:]:
        try:
            nums = [int(s.strip()) for s in row[8]]
        except ValueError:
            Pass = False
            error()
            print("QualityOffset should be integers. Please change this and try again!")
    return Pass


def designStructureCheck(stringArray):
    '''
    Checks if the entries in a design file are {0, 1, 2} and that the first column name is "Sample"
    '''
    Pass = True
    entryOps = [0,1,2]
    if not stringArray[0][0] == "Sample":
        Pass = False
        error()
        print("The name of the first column of the design file should be 'Sample'. Please fix this and try again!")
    for row in stringArray[1:]:
        try:
            nums = [int(s.strip()) for s in row[1:]]
            if not all([n in entryOps for n in nums]):
                Pass = False
                error()
                print(str(nums) + " contains values that are not permitted (" + str(entryOps) + "). Please fix this and try again!")
        except ValueError:
            Pass = False
            error()
            print(str(s) + " is not an integer. Entries in design file matrix should be integers. Please change this and try again!")
    return Pass



def outputPass(checks, filename):
    '''
    Checks if all tests have passed
    '''
    if all(checks.values()):
        print('''
            ---------------------------------------------------------------------------\n
                                    Your file has passed the check!                       
            ---------------------------------------------------------------------------\n
            ''')
    else:
        for key, value in checks.items():
            if not value:
                print(key + " has failed the test. Please fix this before launching the pipelines.")



def readsetValidator(readsetFile):
    '''
    Validator for readset file. activated by -r
    '''
    readsetDict = {}
    readsetStr = loadFile(readsetFile)
    print("\n\nAnalyzing Readset File "  + readsetFile + ":\n")
    readsetDict['unicodePass'] = unicodeCheck(readsetStr)
    readsetDict['DosNewlinePass'] = DosNewline(readsetStr)
    readsetDict['headerPass'] = readsetHeaderCheck(readsetStr[0])
    readsetDict['rowSizePass'] = rowSizeCheck(readsetStr)
    readsetDict['trailingSpacesPass'] = trailingSpacesCheck(readsetStr)
    readsetDict['readsetRunTypePass'] = readsetRunTypeCheck(readsetStr)
    readsetDict['AdapterPass'] = readsetAdapterCheck(readsetStr)
    readsetDict['QualityOffsetPass'] = readsetIntegerCheck(readsetStr)
    print("\n\nSummary Report for file: "  + readsetFile + ":\n")
    outputPass (readsetDict, readsetFile)
    


def designValidator(designFile):
    '''
    Validator for design file. activated by -d
    '''
    designDict = {}
    designStr = loadFile(designFile)
    print("\n\nnalyzing Design File "  + designFile + ":\n")
    designDict['unicodePass'] = unicodeCheck(designStr)
    designDict['DosNewlinePass'] = DosNewline(designStr)
    designDict['rowSizePass'] = rowSizeCheck(designStr)
    designDict['trailingSpacesPass'] = trailingSpacesCheck(designStr)
    designDict['StructurePass'] = designStructureCheck(designStr)
    print("\n\nSummary Report for file: "  + designFile + ":\n")
    outputPass (designDict, designFile)



## main:
## parse command line options:

designFile = readsetFile = None
try:
    options,remainder = getopt.getopt(sys.argv[1:],'r:d:h',["readsetFile=","designFile=","help"])
except getopt.GetoptError, e:
    print("Error - "+str(e)+". See help ('-h' or '--help')")
    sys.exit(2)

for opt,arg in options:
    if(opt in ["-r","--readset"]):
        readsetFile=arg
    elif(opt in ["-d","--design"]):
        designFile=arg
    elif(opt in ["-h","--help"]):
        help()
        sys.exit()


## check if there is a single input file and if it is a valid file:

if (readsetFile == None and designFile == None):
    error()
    print("You need to provide a readset file (-r) or a design file (-d) to be validated!")
    help()
    sys.exit(2)
elif (readsetFile != None):
    if not os.path.isfile(readsetFile):
        error()
        print("The readset file provided does not exist. Please provide a valid readset file!")
        sys.exit(2)
    readsetValidator(readsetFile)
elif (designFile != None):
    if not os.path.isfile(designFile):
        error()
        print("The design file provided does not exist. Please provide a valid design file!")
        sys.exit(2)
    designValidator(designFile)









