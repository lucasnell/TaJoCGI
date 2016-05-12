
cimport cython


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
cpdef winAnalyze_c(list knownWin):
    
    cdef:
        float g, c, cpg, gc, obsExp
        size_t i, n = len(knownWin)
        bytes wini, winii
    
    g = 0
    c = 0
    cpg = 0
    for i in range(n):
        wini = knownWin[i]
        if wini == b'G':
            g += 1
            if i > 0:
                winii = knownWin[(i - 1)]
                if winii == b'C':
                    cpg += 1
            # except IndexError: pass
        elif wini == b'C':
            c += 1
    
    gc = (g + c) / n
    
    if gc > 0:
        obsExp = (cpg / n) / (gc / 2)**2
    else:
        obsExp = 0.0
    
    return gc, obsExp


