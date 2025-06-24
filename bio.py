import pandas
import numpy as np 
import math

UP_ARROW = "\u2191"
RIGHT_ARROW = "\u2192"
DOWN_ARROW = "\u2193"
LEFT_ARROW = "\u2190"
DOWN_RIGHT_ARROW = "\u2198"
UP_LEFT_ARROW = "\u2196"
STOP_ARROW = "\u002e"
GAP_ARROW = "\u002D"

def readSubsMat(fn: str) -> dict:
    smdf = pandas.read_csv(fn, sep='\\s+', comment='#', index_col=0)
    na = smdf.to_numpy()
    smdict = {}
    for k, v in smdf.iterrows():
        ssdict = {}
        for i, j in v.items():
            ssdict[i] = j
        smdict[k] = ssdict
    return smdict

def readFasta(fn: str) -> dict:
    f = open(fn, 'r')
    seqCollection = {}
    acc = ''
    line = f.readline()
    while(line):
        if line.startswith('>'):
            acc = line.split(' ')[0]
            acc = acc[1:].strip()
            seqCollection[acc] = ''
        else:
            if len(acc) > 0:
                seqCollection[acc] = seqCollection[acc] + line.strip()
        line = f.readline()
    return seqCollection

def water(s1: str, s2:str, gapPenalty:int,  sm:dict) -> dict:
    nRows = len(s1) + 1
    nCols = len(s2) + 1

    scoreMatrix = np.full([nRows, nCols], 0)
    traceback = np.full([nRows, nCols], '.')

    maxScore = -1
    maxScoreIndex = (-1,-1)
    for i in range(1, nRows):
        for j in range(1, nCols):
            matchScore = sm[s1[i-1]][s2[j-1]]
            diagonalScore = scoreMatrix[i-1, j-1] + matchScore
            upScore = scoreMatrix[i-1, j] + gapPenalty
            leftScore = scoreMatrix[i, j-1] + gapPenalty

            score = max(0, diagonalScore, upScore, leftScore)
            scoreMatrix[i, j] = score

            if score == 0:
                traceback[i, j] = STOP_ARROW
            elif score == leftScore:
                traceback[i, j] = LEFT_ARROW
            elif score == upScore:
                traceback[i, j] = UP_ARROW
            elif score == diagonalScore:
                traceback[i, j] = UP_LEFT_ARROW

            if score >= maxScore:
                maxScore = score
                maxScoreIndex = (i, j)
    
    s1align = ''
    s2align = ''
    consensus = ''

    i, j = maxScoreIndex
    alignScore = 0

    s1pos = []
    s2pos = []
    while traceback[i, j] != STOP_ARROW:
        s1char = s1[i - 1]
        s2char = s2[j - 1]
        if traceback[i, j] == UP_LEFT_ARROW:
            s1align = s1char + s1align
            s1pos.append(i-1)
            s2align = s2char + s2align
            s2pos.append(j-1)
            if s1char == s2char:
                consensus = '|' + consensus
            else:
                consensus = '.' + consensus
            i = i - 1
            j = j - 1
            alignScore = alignScore + sm[s1char][s2char]
        elif traceback[i, j] == UP_ARROW:
            s1align = s1char + s1align
            s2align = '-' + s2align
            consensus = ' ' + consensus
            i = i - 1
            alignScore = alignScore - 1
        elif traceback[i, j] == LEFT_ARROW:
            s1align = '-' + s1align
            s2align = s2char + s2align
            consensus = ' ' + consensus
            j = j - 1
            alignScore = alignScore - 1

    return {
        's1': s1align,
        's1start': min(s1pos),
        's1end': max(s1pos),
        's2': s2align,
        's2start': min(s2pos),
        's2end': max(s2pos),
        's2': s2align,
        'consensus': consensus,
        'score': alignScore
    }

def needle(s1:str, s2:str, gapPenalty:int, sm:dict) -> dict:
    nRows = len(s1) + 1
    nCols = len(s2) + 1

    scoreMatrix = np.full([nRows, nCols], 0)
    traceback = np.full([nRows, nCols], '-')
    arrow = '.'

    for row in range(nRows):
        for col in range(nCols):
            if row == 0 and col == 0:
                score = 0
                arrow = GAP_ARROW
            elif row == 0:
                prevScore = scoreMatrix[row, col-1]
                score = prevScore + gapPenalty
                arrow = LEFT_ARROW
            elif col == 0:
                prevScore = scoreMatrix[row-1, col]
                score = prevScore + gapPenalty
                arrow = UP_ARROW
            else:
                leftScore =scoreMatrix[row, col-1] + gapPenalty
                upScore = scoreMatrix[row-1, col] + gapPenalty
                
                matchScore = sm[s1[row-1]][s2[col-1]]
                diagonalLeftScore = scoreMatrix[row-1, col-1] + matchScore

                score = max(leftScore, upScore, diagonalLeftScore)

                if score == leftScore:
                    arrow = LEFT_ARROW
                elif score == upScore:
                    arrow = UP_ARROW
                elif score == diagonalLeftScore:
                    arrow = UP_LEFT_ARROW
            
            traceback[row, col] = arrow
            scoreMatrix[row, col] = score

    s1align = ''
    s2align = ''
    consensus = ''

    i = len(s1)
    j = len(s2)
    alignScore = 0

    while traceback[i, j] != GAP_ARROW:
        s1char = s1[i - 1]
        s2char = s2[j - 1]
        if traceback[i, j] == UP_LEFT_ARROW:
            s1align = s1char + s1align
            s2align = s2char + s2align
            if s1char == s2char:
                consensus = '|' + consensus
            else:
                consensus = '.' + consensus
            i = i - 1
            j = j - 1
            alignScore = alignScore + sm[s1char][s2char]
        elif traceback[i, j] == UP_ARROW:
            s2align = '-' + s2align
            s1align = s1char + s1align
            consensus = ' ' + consensus
            i = i - 1
            alignScore = alignScore - 1
        elif traceback[i, j] == LEFT_ARROW:
            s1align = '-' + s1align
            s2align = s2char + s2align
            consensus = ' ' + consensus
            j = j - 1
            alignScore = alignScore - 1

    return {
        's1': s1align,
        's2': s2align,
        'consensus': consensus,
        'score': alignScore
    }

def complement(s:str) -> str:
    c = ''
    iubc = {
        'A': 'T',
        'C': 'G',
        'G': 'C',
        'T': 'A',
        'M': 'K', # aMino -> Keto
        'R': 'Y', # puRine -> pYrimidine
        'W': 'S', # Weak -> Strong
        'S': 'W', # Strong -> Weak
        'Y': 'R', # pYrimidine -> puRine
        'K': 'M', # Keto -> aMino
        'V': 'B', # not T -> not A
        'H': 'D', # not G -> not C
        'D': 'H', # not C -> not G
        'B': 'V', # not A -> not T
        'N': 'N'
    }
    for n in s:
        c = c + iubc[n]
    return c
        
def reverse(s:str) -> str:
    return s[::-1]

def revcomp(s: str) -> str:
    return reverse(complement(s))

def oligonn(oligo:str, conc:float, therm:pandas.DataFrame) -> dict:
    oligo = oligo.upper()
    R = 1.99 # cal/K/mol
    dH = therm.at['init', 'dH']
    dG = therm.at['init', 'dG']
    dS = therm.at['init', 'dS']
    for i in range(0, len(oligo)-1):
        dH = dH + therm.at[oligo[i:i+2], 'dH']
        dG = dG + therm.at[oligo[i:i+2], 'dG']
        dS = dS + therm.at[oligo[i:i+2], 'dS']
    return {
        'dH': float(dH),
        'dG': float(dG),
        'dS': float(dS),
        'Tm': float(dH * 1000 / (dS + (R * math.log(conc/4))) - 273.15)
    }
