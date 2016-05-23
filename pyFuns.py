
from collections import deque


class CGIPassThroughError(LookupError):
    """Raise this to indicating a specific event, depending on the function."""


def winAnalyze(knownWin):
    """GC content and Obs_CpG/Exp_CpG for a single window.
    obsExp formula is from 10.1073/pnas.0510310103.
    """
    
    g = knownWin.count('G')
    c = knownWin.count('C')
    n = len(knownWin)
    gc = (g + c) / n
    cpg = ''.join(knownWin).count('CG') / n
    try:
        obsExp = cpg / (gc / 2) ** 2
    except ZeroDivisionError:
        obsExp = 0.0
    return gc, obsExp



def newWinSearch(fastaList, startPos, winSize):
    
    """Searches for new window that match criteria for TJa; returns start position.

    Args:
        fastaList (list): List derived from a single chromosome's fasta file.
        startPos (int): Position to start looking for matches.
        winSize (int): Length of window.

    Returns:
        int: Starting position of match. If it returns 9000000000000, it hit end of 
            chromosome without matching.
    """
    
    i = startPos + winSize
    win = deque(fastaList[startPos:i])
    gc, obsExp = winAnalyze(list(win))
    while (gc < 0.55) | (obsExp < 0.65):
        try:
            win.append(fastaList[i])
        except IndexError:
            raise CGIPassThroughError('End of chromosome reached.')
        win.popleft()
        gc, obsExp = winAnalyze(list(win))
        i += 1
    return i - winSize




def rollBack(fastaList, startPos, winSize):
    
    """After a non-match, this rolls back by 1 bp until criteria met; returns new end.

    Args:
        fastaList (list): List derived from a single chromosome's fasta file.
        startPos (int): Starting position of the non-match sequence.
        winSize (int): Length of window.

    Returns:
        int: Ending position of match (ending is non-inclusive, so it can be 
            more concisely used for ranges).
    """
    
    s = startPos
    
    win = deque(fastaList[s:(s + winSize)])
    gc, obsExp = winAnalyze(list(win))
    
    i = 1
    
    while ((gc < 0.55) | (obsExp < 0.65)) & (i < winSize):
        s = startPos - i
        win.appendleft(fastaList[s])
        win.pop()
        gc, obsExp = winAnalyze(list(win))
        i += 1
    
    if (gc >= 0.55) & (obsExp >= 0.65):
        return s + winSize
    
    return startPos


def shrinkIsland(fastaList, startPos, endPos, innerMin = 200):
    
    """Shrink a sequence by 1 bp on each side until it meets TJa criteria.

    Args:
        fastaList (list): List derived from a single chromosome's fasta file.
        startPos (int): Start point for sequence to be shrunk. 
        endPos (int): End point for sequence to be shrunk; it should be non-inclusive, 
            so it can be more concisely used for ranges.
        innerMin (int): 'Inner minimum', or the minimum length required inside functions.
            This contrasts the 500 bp set by the TJa, and is included because the TJa
            paper appears to include the 500 bp limit at the very end. Defaults to 200.

    Returns:
        tuple: Two ints, starting position of match, and ending position of match
            (ending is non-inclusive, so it can be more concisely used for ranges).
    """
    
    newStart, newEnd = (startPos, endPos)
    
    gc, obsExp = winAnalyze(fastaList[newStart:newEnd])
    
    while (gc < 0.55) | (obsExp < 0.65):
        newStart += 1
        newEnd -= 1
        if (newEnd - newStart) < innerMin:
            raise CGIPassThroughError("Shrunk below innerMin.", endPos)
        gc, obsExp = winAnalyze(fastaList[newStart:newEnd])
    
    return newStart, newEnd






def oneIsland(fastaList, startPos, lastStart, winSize = 200):
    
    """Search for and, if found, analyze one CpG island.
    
    Args:
        fastaList (list): MV derived from a single chromosome's fasta file.
        startPos (int): Position at which to start looking for the island.
        lastStart (int): Last possible starting position.
        winSize (int): Length of sliding window. Defaults to 200.
    
    Returns:
        tuple: Start and stop points for the CpG island. Use the stop point in the 
            tuple as the next position at which to search for a CpG island.
    
    Special cases:
        1) Raises CGIPassThroughError("Shrunk below innerMin.", slidInd) if the shrinking
           step invalidates the putative CpG island. Don't add this tuple to the list 
           of all CpG islands, but use `slidInd` as new starting point.
        2) Raises CGIPassThroughError('End of chromosome reached.') if it reaches the 
           end of the chromosome w/o finding a CpG island. Break any loops when this 
           occurs.
    """
    
    # Putative new island start
    putStart = newWinSearch(fastaList, startPos, winSize)
    
    # Bc we know they meet criteria, we just quickly initialize them as high values.
    gc = obsExp = 1.0
    
    # Getting start position of next window that no longer meets criteria.
    nonMatchStart = putStart
    while (nonMatchStart < lastStart) & (gc >= 0.55) & (obsExp >= 0.65):
        nonMatchStart += winSize
        # Don't allow to reach beyond end of chromosome.
        nonMatchStart = min([nonMatchStart, lastStart])
        win = fastaList[nonMatchStart:(nonMatchStart + winSize)]
        gc, obsExp = winAnalyze(win)
    
    # Ending position of rolling back by 1 bp until it meets the criteria
    rbEnd = rollBack(fastaList, nonMatchStart, winSize)
    
    # If this new, larger strand doesn't meet criteria, it's shrunk by 1bp per side
    #   until it does.
    shrunkStart, shrunkEnd = shrinkIsland(fastaList, putStart, rbEnd)
    
    return shrunkStart, shrunkEnd






def oneChrom(fastaList, winSize = 200, minSep = 100, minLen = 500):
    
    """Slide through one chromosome looking for CpG islands.
    
    Args:
        fastaList (list): List derived from a single chromosome's fasta file.
        winSize (int): Length of window to slide. Defaults to 200.
        minSep (int): Minimum distance between islands, otherwise we combine. 
            Defaults to 100.
        minLen (int): Minimum length of islands. Defaults to 500.
    
    Returns:
        list: Set of start and stop points for all putative CpG islands.
    """
    
    allCGI = []
    longCGI = []
    
    newIslInd = 0
    
    lastStart = (len(fastaList) - winSize)
    
    while newIslInd <= lastStart:
        try:
            cgiStart, cgiEnd = oneIsland(fastaList, newIslInd, lastStart)
        except CGIPassThroughError as f:
            # See docstring for `oneIsland` for more info on these errors.
            if 'Shrunk below innerMin.' in f.args:
                newIslInd = f.args[1]
                continue
            elif 'End of chromosome reached.' in f.args:
                break
            else:
                raise NotImplementedError('Unknown message for CGIPassThroughError.')
        
        # If it's within `minSep` of previous CGI --> combine + shrink
        # This is inside try-except bc it'll throw NameError if it indexes `allCGI` using
        #   `lastRowCGI`; this is faster way to prevent it doing this check on first CGI
        try:
            prevRow = allCGI[lastRowCGI]
        except NameError:
            pass
        else:
            if (cgiStart - prevRow[1] + 1) <= minSep:
                try:
                    shrunkStart, shrunkEnd = shrinkIsland(fastaList, prevRow[0], cgiEnd)
                # If new, larger CGI doesn't hold up to shrinking, keep old CGI + look
                #   for the next CGI > 100bp from previous
                except CGIPassThroughError:
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




