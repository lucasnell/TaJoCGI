
cimport cython

from cython cimport view
import numpy as np
cimport numpy as np


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
cpdef tuple winAnalyze(list knownWin):
    
    cdef:
        float g, c, cpg, gc, obsExp
        size_t i, n = len(knownWin)
        str wini, winii
    
    g = 0
    c = 0
    cpg = 0
    for i in range(n):
        wini = knownWin[i]
        if wini == 'G':
            g += 1
            if i > 0:
                winii = knownWin[(i - 1)]
                if winii == 'C':
                    cpg += 1
        elif wini == 'C':
            c += 1
    
    gc = (g + c) / n
    
    if gc > 0:
        obsExp = (cpg / n) / (gc / 2)**2
    else:
        obsExp = 0.0
    
    return gc, obsExp



@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
cdef int [:] strList_to_intMV(list str_list):
    """Convert from string list to integer MemoryView."""
    cdef size_t i, n = len(str_list)
    empty_arr = view.array(shape = (n,), itemsize = sizeof(int), format = "i")
    cdef int [:] int_mv = empty_arr
    cdef int i_int
    
    for i in range(n):
        i_int = ord(str_list[i])
        int_mv[i] = i_int
    
    return int_mv

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
cdef int [:] intMV_to_strList(int [:] int_mv):
    """Convert from integer MemoryView to string list."""
    cdef size_t i, n = len(int_mv)
    cdef list str_list = []
    cdef str i_str
    
    for i in range(n):
        i_str = chr(int_mv[i])
        str_list += [i_str]
    
    return str_list

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
cdef tuple winAnalyze_c(int [:] knownWin):
    
    """GC content and Obs_CpG/Exp_CpG for a single window.
    
    obsExp formula is from 10.1073/pnas.0510310103."""
    
    cdef:
        float g = 0, c = 0, cpg = 0, gc, obsExp
        size_t i, n = len(knownWin)
        int wini, winii, g_int = 71, c_int = 67
        float n_float = <float>n
    
    # 'G' = 71, 'C' = 67
    
    for i in range(n):
        wini = knownWin[i]
        if wini == g_int:
            g += 1
            if i > 0:
                winii = knownWin[(i - 1)]
                if winii == c_int:
                    cpg += 1
        elif wini == c_int:
            c += 1
    
    gc = (g + c) / n_float
    
    if gc > 0:
        obsExp = (cpg / n_float) / (gc / 2)**2
    else:
        obsExp = 0.0
    
    return gc, obsExp

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
cdef size_t newWinSearch_c(int [:] fastaMV, size_t startPos, size_t winSize):
    
    """Searches for new window that match criteria for TJa and returns start position.
    
    Args:
        fastaMV (MemoryView): MV derived from a single chromosome's fasta file.
        startPos (size_t): Position to start looking for matches.
        winSize (size_t): Length of window.
    
    Returns:
        (size_t): Starting position of match.
    
    Example:
        s = newWinSearch(chrFasta, 0, 200)
        # To extract the sequence from 'chrFasta'
        win = chrFasta[s:(s+winSize)]
    """
    
    cdef:
        size_t new_start = startPos, new_end = startPos + winSize, n = len(fastaMV)
        int [:] win
        float gc = 0, obsExp = 0
        
    while (gc < 0.55) | (obsExp < 0.65):
        
        if new_end > n:
            return 9000000000000
        
        win = fastaMV[new_start:new_end]
        gc, obsExp = winAnalyze_c(win)
        new_end += 1
        new_start += 1
    
    new_start -= 1
    
    return new_start



@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
cdef size_t rollBack_c(int [:] fastaMV, size_t startPos, size_t winSize):
    
    """After a non-match, this rolls back by 1 bp until criteria met; returns new end.
    
    Args:
        fastaMV (MemoryView): MV derived from a single chromosome's fasta file.
        startPos (size_t): Starting position of the non-match sequence.
        winSize (size_t): Length of window.
    
    Returns:
        (size_t): *Ending* position of match (ending is non-inclusive, so it can be 
            more concisely used for ranges).
    """
    
    cdef:
        size_t new_startPos = startPos, lastStart = startPos - winSize + 1
        int [:] win
        float gc = 0, obsExp = 0
    
    while (gc < 0.55) | (obsExp < 0.65):
        new_startPos -= 1
        # If it goes beyond `lastStart`, return the end position for 
        #   the last matching window.
        if new_startPos < lastStart:
            return startPos
        win = fastaMV[new_startPos:(new_startPos + winSize)]
        gc, obsExp = winAnalyze_c(win)
    
    return new_startPos + winSize


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
cdef tuple shrinkIsland_c(int [:] fastaMV, size_t startPos, size_t endPos, 
                          size_t innerMin = 200):
    
    """Shrink seq. by 1 bp on each side until criteria met; return new start + end.

    Args:
        fastaMV  (MemoryView): MV derived from a single chromosome's fasta file.
        startPos (size_t): Starting position of the non-match sequence.
        endPos (size_t): Ending position of the non-match sequence; should be 
            non-inclusive for making ranges easier.
        innerMin (size_t): 'Inner minimum', or the minimum length required inside 
            functions. This contrasts the 500 bp set by the TJa, and is included 
            because the TJa paper appears to include the 500 bp limit at the very end. 
            Defaults to 200.

    Returns:
        (size_t): Two ints, starting position of match, and ending position of match
            (end is non-inclusive, so it can be more concisely used for ranges).
    """
    
    cdef:
        size_t new_start = startPos, new_end = endPos
        int [:] win
        float gc = 0, obsExp = 0
    
    while (gc < 0.55) | (obsExp < 0.65):
        if (new_end - new_start) < innerMin:
            return 0, 0
        win = fastaMV[new_start:new_end]
        gc, obsExp = winAnalyze_c(win)
        new_start += 1
        new_end -= 1
    
    new_start -= 1
    new_end += 1
    
    return new_start, new_end




@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
cdef tuple oneIsland(int [:] fastaMV, size_t startPos, size_t lastStart, size_t n, 
                     size_t winSize = 200):
    
    """Search for and analyze one CpG island.
    
    Args:
        fastaMV (MemoryView): MV derived from a single chromosome's fasta file.
        startPos (size_t): Position at which to start looking for the island.
        lastStart (size_t): Last possible starting position.
        n (size_t): Length of fastaMV.
        winSize (size_t): Length of sliding window. Defaults to 200.

    Returns:
        tuple: Start and stop points for the CpG island. Use the second item in the 
            tuple to search for the next CpG island.
            
            Special return cases:
            1) Returns (9000000000000, slidInd) if the shrinking step invalidates the 
                putative CpG island. Don't add this tuple to the list of all CpG islands.
            2) Returns (0, 9000000000000) if it reaches the end of the chromosome w/o 
                finding a CpG island.
    """
    cdef:
        size_t nonMatchStart, rbEnd
        int [:] win
        size_t putStart, shrunkStart, shrunkEnd
        float obsExp, gc
    
    # Putative new island start
    putStart = newWinSearch_c(fastaMV, startPos, winSize)
    # `newWinSearch_c` returns 9000000000000 when it reaches the end of the chrom 
    #   w/o finding CGI
    if putStart == 9000000000000:
        return 0, 9000000000000
    
    # Bc we know they meet criteria, we just quickly initialize them as high values.
    gc = obsExp = 1.0
    
    # Getting start position of next window that no longer meets criteria.
    nonMatchStart = putStart
    while (nonMatchStart < lastStart) & (gc >= 0.55) & (obsExp >= 0.65):
        nonMatchStart += winSize
        # Don't allow to reach beyond end of chromosome.
        if nonMatchStart > lastStart:
            nonMatchStart = lastStart
        win = fastaMV[nonMatchStart:(nonMatchStart + winSize)]
        gc, obsExp = winAnalyze_c(win)
    
    # Ending position of rolling back by 1 bp until it meets the criteria
    rbEnd = rollBack_c(fastaMV, nonMatchStart, winSize)
    
    # If this new, larger strand doesn't meet criteria, it's shrunk by 1bp per side
    #   until it does.
    shrunkStart, shrunkEnd = shrinkIsland_c(fastaMV, putStart, rbEnd)
    # `shrinkIsland_c` returns (0,0) when it shrinks to < 200 bp and 
    #   criteria aren't met.
    if shrunkStart == shrunkEnd == 0:
        return 9000000000000, rbEnd
    
    return shrunkStart, shrunkEnd




cpdef list oneChrom(list fastaList, size_t winSize = 200, size_t minSep = 100, 
                    size_t minLen = 500):
    
    """Slide through one chromosome looking for CpG islands.
    
    Args:
        fastaList (list): List derived from a single chromosome's fasta file.
        winSize (size_t): Length of window to slide. Defaults to 200.
        minSep (size_t): Minimum distance between islands, otherwise we combine. 
            Defaults to 100.
        minLen (size_t): Minimum length of islands. Defaults to 500.
    
    Returns:
        list: Set of start and stop points for all putative CpG islands.
    """
    cdef:
        int [:] fastaMV = strList_to_intMV(fastaList)
        list allCGI = [], longCGI = []
        long [:] prevRow = np.zeros(2, int)
        size_t newIslInd = 0, n = len(fastaMV), lastStart = (n - winSize)
        size_t cgiStart, cgiEnd, shrunkStart, shrunkEnd, lastRowCGI, i
    
    while newIslInd <= lastStart:
        cgiStart, cgiEnd = oneIsland(fastaMV, newIslInd, lastStart, n)
        # The following if/else statements account for special returns of `oneIsland`:
        #   - Breaking when it reaches the chromosome end and returns (0, 9000000000000)
        #   - Continuing onto next iteration when shrinking invalidates putatitive CGI
        if cgiEnd == 9000000000000:
            break
        if cgiStart == 9000000000000:
            newIslInd = cgiEnd
            continue
        # If it's within `minSep` of previous CGI --> combine + shrink
        if len(allCGI) > 0:
            prevRow[0] = allCGI[lastRowCGI][0]
            prevRow[1] = allCGI[lastRowCGI][1]
            if (cgiStart - prevRow[1] + 1) <= 100:
                shrunkStart, shrunkEnd = shrinkIsland_c(fastaMV, prevRow[0], 
                                                        cgiEnd)
                # If new, larger CGI doesn't hold up to shrinking, keep old CGI + look for 
                #   the next CGI > 100bp from previous
                if shrunkStart == shrunkEnd == 0:
                    newIslInd += 101
                    continue
                # If it does hold up, remove last entry in `allCGI` before appending
                allCGI = allCGI[:-1]
                cgiStart, cgiEnd = shrunkStart, shrunkEnd
        allCGI += [(cgiStart, cgiEnd)]
        lastRowCGI = len(allCGI) - 1
        newIslInd = cgiEnd
        
    # Now go back through and make sure all CGI are >= `minLen` bp long
    for i in range(len(allCGI)):
        if (allCGI[i][1] - allCGI[i][0]) >= minLen:
            longCGI += [allCGI[i]]
    
    return longCGI


