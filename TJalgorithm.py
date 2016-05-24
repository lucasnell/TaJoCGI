#!/usr/local/apps/anaconda/3-2.2.0/bin/python

"""
Implementation of Takai and Jones’ (2002) algorithm (TJa) for finding CpG islands.
"""


import numpy as np
import re
from multiprocessing import Pool
import gzip
import argparse as ap
import mmap
import cyFuns as cy
# If you want to use full Python version, uncomment below and replace all 'cy.' with 'py.'
# import pyFuns as py




# ==========================================
# ==========================================

#   Functions

# ==========================================
# ==========================================





def readFastaNames(fastaFile, oneSeq = False):
    
    """Read an entire fasta file, and return list of sequence name(s).

    Args:
        fastaFile (str): Name of fasta file. Accepts uncompressed and gzipped files.
        oneSeq (bool): Is there only one sequence present in this file? Defaults to False.

    Returns:
        list: List of sequence names within fasta file.
    """

    assert all([fastaFile.__class__ == str, oneSeq.__class__ == bool]), \
        '`fastaFile` must be a string and `oneSeq` must be boolean.'
    
    seqNames = []
    if fastaFile.endswith('gz'):
        file = gzip.open(fastaFile, 'r')
    else:
        file = open(fastaFile, 'r')

    file_mm = mmap.mmap(file.fileno(), 0, access = mmap.ACCESS_READ)
    
    if oneSeq:
        ind = file_mm.find(b">")
        file_mm.seek(ind)
        seqName = re.sub('(\n|> |>)', '',
                         file_mm.readline().decode()
                         ).split(' ')[0]
        seqNames += [seqName]
    else:
        while True:
            ind = file_mm.find(b">")
            try:
                file_mm.seek(ind)
            except ValueError:
                break
            else:
                seqName = re.sub('(\n|> |>)', '', 
                                 file_mm.readline().decode()
                                 ).split(' ')[0]
                seqNames += [seqName]
            
    file_mm.close()
    file.close()
    
    assert seqNames.__class__ == list, 'Function returned a non-list.'
    assert all([x.__class__ == str for x in seqNames]), \
        'Item(s) in return list are non-strings.'

    return seqNames




def readSeq(fastaFile, seqName):
    
    """Read a single sequence from a fasta file, and return list of sequence.
    
    Args:
        fastaFile (str): Name of fasta file. Accepts uncompressed and gzipped files.
        seqName (str): Name of sequence within fasta file.
    
    Returns:
        list: Sequence split into a list of entirely uppercase letters
            (e.g., 'aCgT' --> ['A', 'C', 'G', 'T']).
    """
    
    assert all([fastaFile.__class__ == str, seqName.__class__ == str]), \
        '`fastaFile` and `seqName` must both be strings.'
    
    if fastaFile.endswith('gz'):
        file = gzip.open(fastaFile, 'r')
    else:
        file = open(fastaFile, 'r')
        
    file_mm = mmap.mmap(file.fileno(), 0, access = mmap.ACCESS_READ)
    
    # Finding the position of the sequence name
    nameInd = file_mm.find(seqName.encode('utf-8'))
    file_mm.seek(nameInd)
    
    # The actual sequence will begin right after the next newline
    start = file_mm.find(b'\n') + 1
    file_mm.seek(start)
    
    # The sequence will end right before the next '>'
    end = file_mm.find(b'>')
    
    # Seek to `start`, read a bytes object of length `end - start`, 
    #   convert to string, remove newlines, and make uppercase
    file_mm.seek(start)
    seqStr = file_mm.read(end - start).decode().replace('\n', '').upper()
    seqList = list(seqStr)
    
    file_mm.close()
    file.close()
    
    assert seqList.__class__ == list, 'Function returned a non-list.'
    assert all([x.__class__ == str for x in seqList]), \
        'Item(s) in return list are non-strings.'

    return seqList






def assembleAllNames(fastaFiles, oneSeq = False):
    
    """Returns a list of (<file>, <sequence>) for all sequences in all input fasta files.
    
    Args:
        fastaFiles (list): All fasta file(s) to process.
    
    Returns:
        (list): Tuple(s) of file and sequence for all sequence(s) in all input fasta 
            file(s).
    """
    
    assert fastaFiles.__class__ == list, '`fastaFiles` must be a list.'
    
    fileSeqList = []
    for f in fastaFiles:
        fSeqs = readFastaNames(f, oneSeq)
        fileSeqList += [(f, x) for x in fSeqs]
    
    assert fileSeqList.__class__ == list, 'Function returned a non-list.'
    assert all([x.__class__ == tuple for x in fileSeqList]), \
        'Item(s) in return list are non-tuples.'
    assert (len(fileSeqList) + len(fastaFiles)) == len(np.unique(fileSeqList)), \
        'Duplicate sequence names.'
    
    return fileSeqList






def getCGI(ID_tup):
    
    """Find CpG islands (CGI) in a sequence within a fasta file.

        Args:
            ID_tup (tuple): Tuple of 
                (1) name of fasta file (accepts uncompressed and gzipped files)
                (2) name of sequence within fasta file

        Returns:
            list: Start and end positions for CpG islands with 0-based indexing.
        """

    assert ID_tup.__class__ == tuple, '`ID_tup` must be tuple'
    
    fastaFile = ID_tup[0]
    seqName = ID_tup[1]
    
    assert all([fastaFile.__class__ == str, seqName.__class__ == str]), \
        '`fastaFile` and `seqName` must both be strings.'
    
    # List derived from sequence in fasta file
    fastaList = readSeq(fastaFile, seqName)
    
    # List of CpG island locations
    CpGislands = cy.oneChrom(fastaList)

    # Adding sequence names to 1st column, to later combine all CGI for all sequences
    wSeqName = [tuple([seqName] + list(x)) for x in CpGislands]
    
    return wSeqName





def main():
    
    """Get CpG island for all fasta sequence(s) and save BED file of locations.
    
    Requires the following objects in global environment:
        - fastaFiles (list)
        - onlySingles (bool)
        - maxCores (int)
        - outFile (str)
    """
    
    assert fastaFiles.__class__ == list, '`fastaFiles` must be list'
    assert onlySingles.__class__ == bool, '`onlySingles` must be boolean'
    assert maxCores.__class__ == int, '`maxCores` must be integer'
    assert outFile.__class__ == str, '`outFile` must be string'
    
    # Names of all sequences from all fasta files:
    allNames = assembleAllNames(fastaFiles, onlySingles)
    
    if maxCores > 1:
        with Pool(processes = maxCores) as pool:
            allSeqList_pool = pool.map(getCGI, allNames)
        allSeqList = [row for seq in allSeqList_pool for row in seq]
    else:
        allSeqList = []
        for f in allNames:
            allSeqList += getCGI(f)
    
    allSeqArray = np.array(allSeqList)
    
    # Save as BED file
    np.savetxt('%s.bed' % outFile, allSeqArray, fmt = '%s',
               delimiter = '\t', comments = '')
    
    return





if __name__ == '__main__':
    
    # ==================
    # Setting up parser
    # ==================
    
    Parser = ap.ArgumentParser(
        description = "Implementation of Takai and Jones' (2002) algorithm for " + \
                      "finding CpG islands.")
    
    Parser.add_argument('-c', '--cores', type = int, default = 1,
                        help = "Maximum number of cores to use.")

    Parser.add_argument('-o', '--outFile', required = True,
                        help = "Name of output BED file, not including '.bed' " + \
                               "extension.")
    
    Parser.add_argument('--onlySingleSeqs', dest = 'onlySingles', action = 'store_true',
                        help = 'Do all input fasta files have only a single ' + \
                               'sequence within? Defaults to False.')
    Parser.set_defaults(onlySingles = False)

    Parser.add_argument('assemblyFiles', nargs = '+',
                        help = "Fasta file(s) to process. They can be uncompressed " + \
                               "or gzipped.")
    
    # ==================
    # Reading the arguments
    # ==================
    args = vars(Parser.parse_args())
    maxCores = args['cores']
    outFile = args['outFile']
    fastaFiles = args['assemblyFiles']
    onlySingles = args['onlySingles']
    if fastaFiles.__class__ == str:
        fastaFiles = [fastaFiles]
    
    
    # ==================
    # Running functions
    # ==================
    main()





