#!/usr/local/apps/anaconda/3-2.2.0/bin/python

"""
Implementation of Takai and Jones’ (2002) algorithm (TJa) for finding CpG islands.
============

Minimums
--------
- Length = 500
- %GC = 55
- (Obs_CpG/ Exp_CpG) = 0.65


Algorithm
--------
1. In 200 bp window at beginning of sequence, get %GC & (Obs_CpG/ Exp_CpG). Shift by 1 bp
    until it meets criteria above.
2. If the window meets the criteria, shift the window 200 bp and then evaluate again.
3. Repeat these 200-bp shifts until the window does not meet the criteria.
4. Shift the last window 1 bp back (toward the 5' end) until it meets the criteria.
5. Evaluate total %GC and Obs_CpG / Exp_CpG for this combined strand.
6. If this large CpG island does not meet the criteria, trim 1 bp from each side until
    it meets the criteria.
7. Two individual CpG islands were connected if they were separated by less than 100 bp.
8. Repeat steps 5-6 to evaluate the new sequence segment until it meets the criteria.
9. Reset start position immediately after the CGI identified at step 8 and go to step 1.

[*From original paper, it seems like 500 bp filtering happened at the end.*]


Example usage
--------

# Using on file `Glazer_chr1.fa` with 12 cores

/path/TJalgorithm.py -c 12 Glazer_chr1.fa 

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
    assert (len(allNames) + len(fastaFiles)) == len(np.unique(fileSeqList)), \
        'Duplicate sequence names.'
    
    return fileSeqList






def getCpGisland(fastaFile, seqName):
    
    """Find CpG islands in a sequence within a fasta file.

        Args:
            fastaFile (str): Name of fasta file. Accepts uncompressed and gzipped files.
            seqName (str): Name of sequence within fasta file.

        Returns:
            list: Start and end positions for CpG islands with 0-based indexing.
        """
    
    assert all([fastaFile.__class__ == str, seqName.__class__ == str]), \
        '`fastaFile` and `seqName` must both be strings.'

    fastaList = readSeq(fastaFile, seqName)
    
    CpGislands = cy.oneChrom(fastaList)

    # Adding sequence names to later combine all CpG island for all sequences
    wSeqName = [tuple([seqName] + list(x)) for x in CpGislands]
    
    return wSeqName

















# npOut = np.array(LISTNAME, dtype = [('seq', object), ('start', int), ('end', int)])






 
fastaFiles = ['chrI.fa', 'chrXIX.fa', 'Glazer_unmasked.fa']

fastaFiles = ['/Volumes/MW_18TB/Lucas_Nell/lan/stick/genome/assem/unmasked/' + x for x
              in fastaFiles]


allNames = assembleAllNames(fastaFiles, False)







scafNames, scafSeqs = readFasta(
    '/Volumes/MW_18TB/Lucas_Nell/lan/stick/genome/assem/Ensembl/rawFasta/scaffolds.fa',
    namesOnly = False)

# min([len(z) for z in scafSeqs.values()])




np.array(scafSeqs[(len(scafSeqs) - 1)])














if __name__ == '__main__':
    
    # ==================
    # Setting up parser
    # ==================
    ScriptDescript = 'Implementation of Takai and Jones’ (2002) algorithm (TJa) for ' + \
                     'finding CpG islands.'
    
    Parser = ap.ArgumentParser(description=ScriptDescript)
    
    Parser.add_argument('-c', '--cores', type = int, default = 1,
                        help = "Maximum number of cores to use.")
    
    Parser.add_argument('assemblyFiles', nargs = '+',
                        help = "Fasta file(s) to process. They can be uncompressed " + \
                               "or gzipped.")
    
    Parser.add_argument('--onlySingleSeqs', dest = 'onlySingles', action = 'store_true',
                        help = 'Do all input fasta files have only a single sequence ' + \
                               'within? Defaults to False.')
    Parser.set_defaults(onlySingles = False)
    
    # ==================
    # Reading the arguments
    # ==================
    args = vars(Parser.parse_args())
    maxCores = args['cores']
    fastaFiles = args['assemblyFiles']
    onlySingles = args['onlySingles']
    if fastaFiles.__class__ == str:
        fastaFiles = [fastaFiles]
    
    
    
    
    
    
    
    
    # # ==================
    # # Running code on arguments
    # # ==================
    # 
    # 
    # if maxCores > 1:
    #     with Pool(processes = maxCores) as pool:
    #         allChromoList = pool.map(getFastaIsls, fastaFullPaths)
    #     allChromoListSqueezed = [row for chromo in allChromoList for row in chromo]
    #     allChrArray = np.array(allChromoListSqueezed)
    # else:
    #     allChromoList = []
    #     for f in fastaFullPaths:
    #         fChrList = getFastaIsls(f)
    #         allChromoList += fChrList
    #     allChrArray = np.array(allChromoList)
    # np.savetxt('CGI_%s' % whichAssembly, allChrArray, fmt = '%s',
    #            delimiter = '\t', header = 'chr\tstart\tend', comments = '')





