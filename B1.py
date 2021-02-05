## Vanda Azevedo 2018

import numpy as np
import requests
import sys
import operator
from MatrixInfo import blosum62 as blosum
from MatrixInfo import pam120 as pam

max_iter=5
sys.setrecursionlimit(2000)

class Align(object):
    def __init__(self, A, B, score):
        self.A=A
        self.B=B
        self.score=score

class Scores(object):
    def __init__(self):
        self.scores=list()
        self.l=0
        self.paths_nr=0

    def add(self, A, B):
        id_score=idScore(A,B)
        self.scores.append(Align(A,B,id_score))
        self.paths_nr+=1

    def printScores(self):
        maxi=max_iter
        self.scores=sorted(self.scores, key=operator.attrgetter('score'), reverse=True)
        if(len(self.scores)<maxi):
            maxi=len(self.scores)
        for i in range(0, maxi):
            finalize(self.scores[i].A, self.scores[i].B)
    
    def clear(self):
        self.scores = list()
        self.l=0
        self.paths_nr=0
    
    def local(self, b):
        self.l=b
    
    def test(self):
        for i in range(0, len(self.scores)):
            print(self.scores[i].score, " ", end='')


stype="protein"
gap = -7  # worse than mismatch
h = -15 # worse than gap
chunk=75
sep = "="*90
mini_sep= "="*50
scores = Scores()
local=0
max_paths=1000
matrix_type="pam"

##### ALGORITHMS

def globalAlign(A, B):
    scores.local(0)
    m = len(A)
    n = len(B)
    min_num=-sys.maxsize
    diag = np.zeros((m+1, n+1))
    top= np.zeros((m+1, n+1))
    left= np.zeros((m+1, n+1))

    for i in range(0, m+1):
        top[i][0]=h+gap*i
        diag[i][0]=min_num
        left[i][0]=min_num
    for j in range(0, n+1):
        left[0][j]=h+gap*j
        diag[0][j]=min_num
        top[0][j]=min_num

    diag[0][0]=0
    top[0][0]=h
    left[0][0]=h

    for i in range(1, m+1):
        for j in range(1, n+1):
            match = diag[i-1][j-1]+S(A[i-1], B[j-1])
            insertTop = top[i-1][j-1]+S(A[i-1], B[j-1])
            insertLeft = left[i-1][j-1]+S(A[i-1], B[j-1])
            diag[i][j] = max(match, insertTop, insertLeft)

            openGapTop = diag[i-1][j]+h+gap
            extendTop = top[i-1][j]+gap
            top[i][j] = max(openGapTop, extendTop)

            openGapLeft = diag[i][j-1]+h+gap
            extendLeft = left[i][j-1]+gap
            left[i][j] = max(openGapLeft, extendLeft)

    getMaxScore(diag, top, left, A, B, "", "", len(A), len(B), 1)
    
    #printScoreMatrix(diag,top,left,A,B)
    scores.printScores()
    #scores.test()
    scores.clear()

def localAlign(A, B):
    scores.local(1)
    m = len(A)
    n = len(B)
    diag = np.zeros((m+1, n+1))
    top= np.zeros((m+1, n+1))
    left= np.zeros((m+1, n+1))

    for i in range(0, m+1):
        diag[i][0] = 0
        top[i][0]=-sys.maxsize
        left[i][0]=-sys.maxsize
    for j in range(0, n+1):
        diag[0][j] = 0
        left[0][j]=-sys.maxsize
        top[0][j]=-sys.maxsize
    maxi=0

    for i in range(1, m+1):
        for j in range(1, n+1):
            matchS = diag[i-1][j-1]+S(A[i-1], B[j-1])
            insertTop = top[i-1][j-1]+S(A[i-1], B[j-1])
            insertLeft = left[i-1][j-1]+S(A[i-1], B[j-1])
            diag[i][j] = max(0, matchS, insertTop, insertLeft)

            openGapTop = diag[i-1][j]+h+gap
            extendTop = top[i-1][j]+gap
            top[i][j] = max(openGapTop, extendTop)

            openGapLeft = diag[i][j-1]+h+gap
            extendLeft = left[i][j-1]+gap
            left[i][j] = max(openGapLeft, extendLeft)
            if(matchS>maxi):
                maxi=matchS

    for i in range(1, m+1):
        for j in range(1, n+1):
            matchS = diag[i-1][j-1]+S(A[i-1], B[j-1])
            if(matchS==maxi):
                getMaxScore(diag, top, left, A, B, "", "", i, j, 0)        
    #printScoreMatrix(diag,top,left,A,B)
    scores.printScores()
    #scores.test()
    scores.clear()

# EMBOSS Needle reads two input sequences and writes their optimal global sequence alignment to file. It uses the Needleman-Wunsch alignment algorithm to find the optimum alignment (including gaps) of two sequences along their entire length.
def needle(A, B, matrix, gapopen, gapext):
    scores.local(0)
    runUrl="https://www.ebi.ac.uk/Tools/services/rest/emboss_needle/run"
    params={'email':"vfba@gmail.com", 
            'title':"needle", 
            'matrix':str(matrix), #EBLOSUM, EPAM
            'gapopen':str(gapopen), #1,5,10,15,20,25,50,100
            'gapext':str(gapext), #0.0005,0.001,0.05,0.1,0.2,0.5,0.6,0.8,1.0,5.0,10.0
            'endweight':'false',
            'format':'pair', #pair,markx0,markx1,markx2,markx3,markx10,srspair,score
            'stype':stype,
            'asequence':A,
            'bsequence':B}
    post=requests.post(url=runUrl, data=params)
    job_id=post.text
    statusUrl="http://www.ebi.ac.uk/Tools/services/rest/emboss_needle/status/"+job_id.replace("\n","")
    get= requests.get(statusUrl)
    print(get.content)
    while(get.text!="FINISHED"):
        get= requests.get(statusUrl)
    embossUrl = "http://www.ebi.ac.uk/Tools/services/rest/emboss_needle/result/" + job_id.replace("\n","") + "/aln"
    get= requests.get(embossUrl)
    return get.text

# EMBOSS Water uses the Smith-Waterman algorithm (modified for speed enhancments) to calculate the local alignment of a sequence to one or more other sequences.
def water(A, B, matrix, gapopen, gapext):
    scores.local(1)
    runUrl="https://www.ebi.ac.uk/Tools/services/rest/emboss_water/run"
    params={'email':"vfba@gmail.com", 
            'title':"water", 
            'matrix':str(matrix), #EBLOSUM, EPAM
            'gapopen':str(gapopen), #1,5,10,15,20,25,50,100
            'gapext':str(gapext), #0.0005,0.001,0.05,0.1,0.2,0.5,0.6,0.8,1.0,5.0,10.0
            'format':'pair', #pair,markx0,markx1,markx2,markx3,markx10,srspair,score
            'stype':stype,
            'asequence':A,
            'bsequence':B}
    post=requests.post(url=runUrl, data=params)
    job_id=post.text
    statusUrl="https://www.ebi.ac.uk/Tools/services/rest/emboss_water/status/"+job_id.replace("\n","")
    get= requests.get(statusUrl)
    print(get.content)
    while(get.text!="FINISHED"):
        get= requests.get(statusUrl)
    embossUrl = "https://www.ebi.ac.uk/Tools/services/rest/emboss_water/result/" + job_id.replace("\n","") + "/aln"
    get= requests.get(embossUrl)
    return get.text


##### SCORING
def getMaxScore(diag, top, left, A, B, AlignA, AlignB, i, j, glob):
    matchS=diag[i][j]
    insertTop=top[i][j]
    insertLeft=left[i][j]
    current = max(matchS, insertTop, insertLeft)
    if(scores.paths_nr>max_paths):
        return

    if(matchS == 0 and glob==0):
        scores.add(AlignA[::-1], AlignB[::-1])
        return

    if(i > 0 and j > 0 and ((current != 0 and glob==0) or glob==1)):
        if current == matchS:
            #print("(",i,",",j,")=match")
            tempA = AlignA + A[i-1]
            tempB = AlignB + B[j-1]
            getMaxScore(diag, top, left, A, B, tempA, tempB, i-1, j-1, glob)
        if current == insertTop:
            #print("(",i,",",j,")=left")
            tempA = AlignA + A[i-1]
            tempB = AlignB + '-'
            getMaxScore(diag, top, left, A, B, tempA, tempB, i-1, j, glob)
        if current == insertLeft:
            #print("(",i,",",j,")=top")
            tempA = AlignA + '-'
            tempB = AlignB + B[j-1]
            getMaxScore(diag, top, left, A, B, tempA, tempB, i, j-1, glob)
            
    if(glob==1 and (i==0 or j==0)):
        while i > 0:
            AlignA += A[i-1]
            AlignB += '-'
            i -= 1
        while j > 0:
            AlignA += '-'
            AlignB += B[j-1]
            j -= 1
        scores.add(AlignA[::-1], AlignB[::-1])
        return

def S(a, b):
    s = a+b
    dic = {'AA': 10, 'AG': -1, 'AC': -3, 'AT': -4,
                'GA': -1, 'GG': 7, 'GC': -5, 'GT': -3,
                'CA': -3, 'CG': -5, 'CC': 9, 'CT': 1,
                'TA': -4, 'TG': -3, 'TC': 1, 'TT': 8}
    if(stype=="dna" and s in dic):
        return dic[s]
    elif(stype=="protein" and matrix_type=="blosum" and (a,b) in blosum):
        return blosum[(a,b)]
    elif(stype=="protein" and matrix_type=="blosum" and (b,a) in blosum):
        return blosum[(b,a)]
    elif(stype=="protein" and matrix_type=="pam" and (a,b) in pam):
        return pam[(a,b)]
    elif(stype=="protein" and matrix_type=="pam" and (b,a) in pam):
        return pam[(b,a)]
    else:
        return gap+1


def printScoreMatrix(diag, top, left, A, B):
    m = len(A)
    n = len(B)
    print("* Matrix DIAGONAL *")
    for i in range(0, m+1):
        for j in range(0, n+1):
            print(diag[i][j], ' ', end='')
        print()
    print()
    print("* Matrix TOP *")
    for i in range(0, m+1):
        for j in range(0, n+1):
            print(top[i][j], ' ', end='')
        print()
    print()
    print("* Matrix LEFT *")
    for i in range(0, m+1):
        for j in range(0, n+1):
            print(left[i][j], ' ', end='')
        print()
    print()


def finalize(AlignA, AlignB):
    # Calculate identity, diag and aligned sequences
    fill = ''
    score = 0
    identity = 0
    gapCount=0
    length=len(AlignA)

    for i in range(0, length):
        # Bases =
        if AlignA[i] == AlignB[i]:
            fill += '|'
            identity = identity + 1
            score += S(AlignA[i], AlignB[i])
        # Bases !=
        elif AlignA[i] != AlignB[i] and AlignA[i] != '-' and AlignB[i] != '-':
            score += S(AlignA[i], AlignB[i])
            fill += ':'
        # Uma e gap
        elif AlignA[i] == '-' or AlignB[i] == '-':
            fill += ' '
            score += gap
            gapCount+=1

    identityP = float(identity) / length * 100
    gapCountP = float(gapCount) / length * 100
    print()
    print("#",mini_sep)
    print("#")
    if(scores.l==1):
        print("# [ LOCAL ]")
    else:
        print("# [ GLOBAL ]")
    print("# Length = ",length)
    print('# Identity = ', identity, "/", length, " (%3.1f" % identityP,"%)")
    print('# Gaps = ', gapCount, "/", length, " (%3.1f" % gapCountP,"%)")
    print('# Score = ', score)
    print("#")
    print("#",mini_sep)
    print()
    i=1
    chunki=75
    prev_chunk=0
    while(prev_chunk<length):
        chunk=chunki*i
        if(chunk>length):
            print(prev_chunk+1,"-",length)
        else:
            print(prev_chunk+1,"-",chunk)
        print(AlignA[prev_chunk:chunk])
        print(fill[prev_chunk:chunk])
        print(AlignB[prev_chunk:chunk])
        print()
        prev_chunk=chunk
        i+=1

def idScore(AlignA, AlignB):
    # Calculate identity, score and aligned sequences
    identity = 0
    length=len(AlignA)

    for i in range(0, length):
        # Bases =
        if AlignA[i] == AlignB[i]:
            identity = identity + 1
    if(length!=0):
        identityP = float(identity) / length * 100
    else: 
        return 0
    return identityP


##### UTILS

def parse(A):
    i = A.index('\n')
    B=A[i:]
    C=B.replace("\n","")
    return C

def extractGene(gene):
    url='http://www.ebi.ac.uk/ena/data/view/%s&display=fasta' %(gene)
    r = requests.get(url, auth=('user', 'pass'))
    seq=str(r.content)
    seq=cleanData(seq)
    return seq

# https://www.uniprot.org/help/accession_numbers
def extractProtein(protein):
    url= "https://www.ebi.ac.uk/proteins/api/proteins/%s" %(protein)
    r = requests.get(url, headers={ "Accept" : "text/x-fasta"})
    seq = str(r.content)
    seq=cleanData(seq)
    return seq

def cleanData(seq):
    seq=seq.replace("\\n","\n")
    seq=seq.replace("b'","")
    seq=seq.replace("'","")
    return str(seq)


##### MAIN

if __name__ == '__main__':

    print("\nSelect one of the options:\n")
    print("1) Default GENES: AAO27473 and OWJ99766")
    print("2) Default PROTEINS: P22692 and P24594")
    print("3) Insert GENE accession codes")
    print("4) Insert PROTEIN accession codes")
    print("5) Insert GENE sequence")
    print("6) Insert PROTEIN sequence")
    print("7) EXIT")
    op = int(input())
    print()
    first=0
    while(op!=7):
        if(first!=0):
            ########### MUDAR AQUI
            print("\nSelect one of the options:\n")
            print("1) Default GENES: AAO27473 and OWJ99766")
            print("2) Default PROTEINS: P22692 and P24594")
            print("3) Insert GENE accession codes")
            print("4) Insert PROTEIN accession codes")
            print("5) Insert GENE sequence")
            print("6) Insert PROTEIN sequence")
            print("7) EXIT")
            op = int(input())
            if(op==7):
                break
            print()
        else:
            first=1

        if(op == 1 or op==3 or op==5):
            stype="dna"
            if(op == 1):
                A="AAO27473" # Homo sapiens (human) partial HERC2
                B="OWJ99766" # Mus musculus Herc2 (Herc2) mRNA, complete cds.:
            elif(op==3):
                A = str(input("Gene 1:\n"))
                B = str(input("Gene 2:\n"))
            else:
                A = str(input("Sequence 1:\n"))
                B = str(input("Sequence 2:\n"))
        else:
            stype="protein"
            if(op == 2):
                A="P22692" # Homo sapiens (human) insulin
                B="P24594" # Rat insulin
            elif(op==4):
                A = str(input("Protein Accession 1:\n"))
                B = str(input("Protein Accession 2:\n"))
            else:
                A = str(input("Sequence 1:\n"))
                B = str(input("Sequence 2:\n"))
        if(op==1 or op==3):
            A=parse(extractGene(A))
            B=parse(extractGene(B))
        elif(op==2 or op==4):
            A=parse(extractProtein(A))
            B=parse(extractProtein(B))
        else:
            pass
        
        #if(op==1 or op==2 or op==3 or op==4):
        print("Insert the following parameters INLINE:")
        print(" matrix [EBLOSUM$, EPAM$]")
        print(" gapopen [1, 5, 10, 15, 20, 25, 50, 100]")
        print(" gapextend [0.0005, 0.001, 0.05, 0.1, 0.2, 0.5, 0.6, 0.8, 1.0, 5.0, 10.0]")
        print()
        matrix, gapopen, gapext=input().split()
        gap=-int(gapext)
        h=-int(gapopen)
        print()
        if "BLOSUM" in matrix:
            matrix_type="blosum"
        print(sep)
        print("\n [GLOBAL] Needleman-Wunsch\n")
        print(sep)
        print()
        print(needle(A,B,matrix,gapopen,gapext))
        print()
        print(sep)
        print("\n [LOCAL] Waterman-Eggert\n")
        print(sep)
        print()
        print(water(A,B,matrix,gapopen,gapext))
        print()
        print(sep)
        print("\n [GLOBAL] Algorithm\n")
        print(sep)
        print()
        globalAlign(A, B)
        print()
        print(sep)
        print("\n [LOCAL] Algorithm\n")
        print(sep)
        print()
        localAlign(A,B)
        print()
