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
from timeit import default_timer as timer
from multiprocessing import Pool
import gzip
import argparse as ap
import cyFuns as cy
import pyFuns as py





if __name__ == '__main__':
    
    # Setting up parser
    ScriptDescript = 'Implementation of Takai and Jones’ (2002) algorithm (TJa) for ' + \
                     'finding CpG islands.'
    
    Parser = ap.ArgumentParser(description = ScriptDescript)
    
    Parser.add_argument('-c', '--cores', type = int, metavar = 'C', required = False,
                        default = 1, help = "Maximum number of cores to use.")
    
    Parser.add_argument('assemblyFiles', metavar = 'aF', nargs = '+',
                        help = "Fasta file(s) to process. They can be uncompressed or " + \
                               "gzipped.")
    
    # Now reading the arguments
    args = vars(Parser.parse_args())
    maxCores = args['cores']
    fastaFiles = args['assemblyFiles']
    if fastaFiles.__class__ == str:
        fastaFiles = [fastaFiles]



# ==========================================
# ==========================================

#   Functions

# ==========================================
# ==========================================



def fastaToList(fastaFile, namesOnly = False):
    
    """Reads a fasta file to a list, and returns chromosome name.

    Note: Makes all bases uppercase."""
    
    chrSeqs = []
    chrNames = []
    if fastaFile.endswith('gz'):
        file = gzip.open(fastaFile, 'rt')
    else:
        file = open(fastaFile, 'rt')
    
    if namesOnly:
        for line in file:
            if line.startswith('>'):
                chrNames += [re.sub('[\n>]', '', line.replace('> ', '>').split(' ')[0])]
        file.close()
        return chrNames
    else:
        for line in file:
            if line.startswith('>'):
                chrNames += [re.sub('[\n>]', '', line.replace('> ', '>').split(' ')[0])]
                try:
                    chrSeqs += [seq]
                except NameError:
                    seq = []
            else:
                seq += list(line.strip().upper())
        chrSeqs += [seq]
    
    file.close()
    
    return chrNames, chrSeqs





def newWinSearch(fastaList, startPos, window):
    
    """Searches for new window that match criteria for TJa.

    Args:
        fastaList (list): List derived from a single chromosome's fasta file.
        startPos (int): Position to start looking for matches.
        window (int): Length of window.

    Returns:
        list: Two ints, starting position of match, and ending position
            (ending is plus 1, so it can be more concisely used for ranges).

    Example:
        s,e = newWinSearch(chrFasta, 0, 200)
        # To extract the sequence from 'chrFasta'
        win = chrFasta[s:e]
    """
    
    i = startPos + window
    win = deque(fastaList[startPos:i])
    gc, obsExp = winAnalyze(list(win))
    while (gc < 0.55) | (obsExp < 0.65):
        try:
            win.append(fastaList[i])
        except IndexError:
            return
        win.popleft()
        gc, obsExp = winAnalyze(list(win))
        i += 1
    return [(i - window), i]



def rollBack(fastaList, startPos, window):
    
    """After finding a non-match, this rolls back by 1 bp until it meets criteria.

    Args:
        fastaList (list): List derived from a single chromosome's fasta file.
        startPos (int): Starting position of the non-match sequence.
        window (int): Length of window.

    Returns:
        list: Two ints, starting position of match, and ending position of match
            (ending is plus 1, so it can be more concisely used for ranges).
    """
    
    s = startPos
    
    win = deque(fastaList[s:(s + window)])
    gc, obsExp = winAnalyze(list(win))
    
    i = 1
    
    while ((gc < 0.55) | (obsExp < 0.65)) & (i < window):
        s = startPos - i
        win.appendleft(fastaList[s])
        win.pop()
        gc, obsExp = winAnalyze(list(win))
        i += 1
    
    if (gc >= 0.55) & (obsExp >= 0.65):
        return [s, (s + window)]
    return





def shrinkIsland(startEnd, fastaList, innerMin = 200):
    
    """Shrink a sequence by 1 bp on each side until it meets TJa criteria.

    Args:
        startEnd (list): Start and end points for sequence to be shrunk (ending
            should be plus 1, so it can be more concisely used for ranges).
        fastaList (list): List derived from a single chromosome's fasta file.
        innerMin (int): 'Inner minimum', or the minimum length required inside functions.
            This contrasts the 500 bp set by the TJa, and is included because the TJa
            paper appears to include the 500 bp limit at the very end. Defaults to 200.

    Returns:
        list: Two ints, starting position of match, and ending position of match
            (ending is plus 1, so it can be more concisely used for ranges).
    """
    
    se = list(startEnd.copy())
    
    gc, obsExp = winAnalyze(fastaList[se[0]:se[1]])
    
    while (gc < 0.55) | (obsExp < 0.65):
        se[0] += 1
        se[1] -= 1
        if (se[1] - se[0]) < innerMin:
            return
        gc, obsExp = winAnalyze(fastaList[se[0]:se[1]])
    
    return se


def getOneChrIsl(fastaList, chrName, window = 200, w = 100):
    
    """Slide through one chromosome looking for CpG islands.

    Note: Takes ~3-4 minutes.

    Args:
        fastaList (list): List derived from a single chromosome's fasta file.
        chrName (str): Name of chromosome.
        window (int): Length of window to slide. Defaults to 200.
        w (int): Minimum distance between islands, otherwise we combine. Defaults to 100.

    Returns:
        list: Set of start and stop points for all putative CpG islands.
    """
    
    # List of all start and end points
    allStartsEnds = []
    
    print('--- Started', chrName)
    
    t0 = timer()
    
    ind = 0
    while ind < len(fastaList):
        # Putative new island range
        putIslandRange = newWinSearch(fastaList, ind, window)
        # Will throw TypeError when it reaches end of chromosome
        # noMatchStart is for looking through new 200 bp blocks until criteria not met
        try: noMatchStart = putIslandRange[1]
        except TypeError: break
        
        gc, obsExp = (1.0, 1.0)
        while (gc >= 0.55) & (obsExp >= 0.65):
            # ZeroDivisionError here means criteria met until end of chromosome
            try:
                gc, obsExp = winAnalyze(fastaList[noMatchStart:(noMatchStart + window)])
                noMatchStart += window
                continue
            except ZeroDivisionError:
                noMatchStart = len(fastaList)
                break
        noMatchStart -= window
        
        # Rolling back by 1 bp until it meets the criteria
        rb = rollBack(fastaList, noMatchStart, window)
        
        # If rb is None (rb[1] returns TypeError), none met the criteria
        # In that case, you want noMatchStart (not noMatchStart-1 bc ends are +1)
        try:
            putIslandRange[1] = rb[1]
        except TypeError:
            putIslandRange[1] = noMatchStart
        
        # Evaluating for new, larger strand
        newRange = shrinkIsland(putIslandRange, fastaList)
        
        # If it's within 'w' bp of last end position, combine them and shrink again
        # IndexError for if it's first island, TypeError for if newRange is None
        try:
            lastStartEnd = allStartsEnds[(len(allStartsEnds)-1)]
            if (newRange[0] - lastStartEnd[1] - 1) < w:
                newRange = shrinkIsland([lastStartEnd[0], newRange[1]], fastaList)
                if newRange is not None:
                    allStartsEnds = allStartsEnds[:-1]
        except (IndexError, TypeError): pass
        
        # newRange will be None if (end - noMatchStart) < innerMin
        try: ind = newRange[1]
        except TypeError: ind = putIslandRange[1]; continue
        
        newRange[1] -= 1
        allStartsEnds += [[chrName] + newRange]
    
    allStartsEnds = [x for x in allStartsEnds if (x[2] - x[1] + 1) >= 500]
    
    print('>>>>>> Finished %s in %.2f seconds' % (chrName, timer() - t0))
    
    return allStartsEnds




def getFastaIsls(fastaFile, window = 200, w = 100):
    
    """Find CpG islands for all chromosomes in a fasta.

    Note: Takes ~3-4 minutes per chromosome.

    Args:
        fastaFile (str): Fasta file in which to find CpG islands.
        window (int): Length of window to slide. Defaults to 200.
        w (int): Minimum distance between islands, otherwise we combine. Defaults to 100.

    Returns:
        list: Set of start and stop points for all putative CpG islands.
    """
    
    chrNames, fastaLists = fastaToList(fastaFile)
    
    assert len(chrNames) == len(fastaLists), \
        "Chromosome names and sequences from fasta file are not the same length."
    
    islList = []
    
    for i in range(len(chrNames)):
        fastaList = fastaLists[i]
        chrName = chrNames[i]
        iChrList = getOneChrIsl(fastaList, chrName, window, w)
        islList += iChrList
    
    return islList





# def findAndSave(fastaPath):
#
#     """Find CpG islands on a chromosome and save table to files."""
#
#     CGI = findIslands(fastaPath)
#
#     np.savetxt('chr%i_CGI.txt' % fastaFile, CGI, fmt = '%i',
#                delimiter = '\t', header = 'start\tend', comments = '')
#     print('Finished saving chromosome %i output' % chromo)
#     return
#
#
#
#
#
#
# def getChrIsl(chromo, cluster = onCluster):
#
#     """Get CpG island (CGI) location file for a single chromosome & add chr number."""
#
#     baseName = './files/IslLocations/indiv_chr/chr%i_CGI.txt'
#     if cluster: baseName = '/lustre1/lan/stick/CGI/indiv_chr/chr%i_CGI.txt'
#     x = np.genfromtxt(baseName % chromo, dtype = 'int',
#                       delimiter = '\t', skip_header = 1)
#     y = np.repeat([chromo], len(x))
#     out = np.c_[ y, x ].copy()
#     return out
#
#
# def combineIslTables(cluster = onCluster):
#
#     """Combine all single-chromosome CpG island files (both strands) into one file."""
#
#     allIsl = np.empty((0,3))
#
#     for c in range(1, 21+1):
#         a = getChrIsl(c)
#         allIsl = np.concatenate([allIsl, a], axis = 0).copy()
#
#     fullPath = './files/IslLocations/CpGislands.txt.gz'
#     if cluster: fullPath = '/lustre1/lan/stick/CGI/CpGislands.txt.gz'
#
#     np.savetxt(fullPath, allIsl, fmt = '%i',
#                delimiter = '\t', header = 'chr\tstart\tend', comments = '')
#
#     print('Finished writing.')
#
#     return


if __name__ ==  '__main__':
    if maxCores > 1:
        with Pool(processes = maxCores) as pool:
            allChromoList = pool.map(getFastaIsls, fastaFullPaths)
        allChromoListSqueezed = [row for chromo in allChromoList for row in chromo]
        allChrArray = np.array(allChromoListSqueezed)
    else:
        allChromoList = []
        for f in fastaFullPaths:
            fChrList = getFastaIsls(f)
            allChromoList += fChrList
        allChrArray = np.array(allChromoList)
    np.savetxt('CGI_%s' % whichAssembly, allChrArray, fmt = '%s',
               delimiter = '\t', header = 'chr\tstart\tend', comments = '')





