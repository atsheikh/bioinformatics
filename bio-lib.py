from collections import Counter
import collections
import numpy as np
import random
from math import log
from sys import argv
import re
from itertools import product 
import copy
import time
import operator


"""
Following code is part of WEEK-1 of 
course BIO-INFORMATICS at COURSERA
"""

def remove_duplicates(values):
    output = []
    seen = set()
    for value in values:
        # If value has not been encountered yet,
        # ... add it to both list and set.
        if value not in seen:
            output.append(value)
            seen.add(value)
    return output

def Text (i, k, DNA):	
	return DNA[i:i+k]

def ReverseComplement(DNA):	
	rc = []
	for i in range (0, len(DNA)):		
		if DNA[i] == "A":
			rc.append("T")
		elif DNA[i] == "T":
			rc.append("A")
		elif DNA[i] == "C":
			rc.append("G")
		elif DNA[i] == "G":
			rc.append("C")
			
	rc = ''.join(rc[::-1])
	
	return rc
		
def PatternCount (DNA, PATTERN):
	count = 0
		
	for i in range (0, len(DNA)-len(PATTERN)+1):		
		if  Text(i, len(PATTERN), DNA) == PATTERN:
			count+=1			
			
	return count

def PatterMatch (DNA, PATTERN):
	index = []
		
	for i in range (0, len(DNA)-len(PATTERN)):		
		if  Text(i, len(PATTERN), DNA) == PATTERN:
			index.append(i)			
			
	return index	

def SymbolToNumber (symbol):
	if symbol == "A":
		return 0
	elif symbol == "C":
		return 1
	elif symbol == "G":
		return 2
	elif symbol == "T":
		return 3

def PatternToNumber (PATTERN):
	if not PATTERN:
		return 0
	
	symbol = PATTERN[-1];
	prefix = PATTERN[0:len(PATTERN)-1]
	
	return 4*PatternToNumber(prefix) + SymbolToNumber(symbol)		

def NumberToSymbol (number):
	if number == 0:		
		return "A"
	elif number == 1:		
		return "C"
	elif number == 2:		
		return "G"
	elif number == 3:		
		return "T"	

def NumberToPattern (index, k):
	
	if k == 1:
		return NumberToSymbol(index)		
	
	prefixIndex = int(index / 4)
	r = int(index % 4)
		
	symbol = NumberToSymbol(r)	
	prefixPattern = NumberToPattern(prefixIndex, k-1)
	
	return (prefixPattern+symbol)
		
def FrequentWords (DNA, k):
	FrequentPatterns = []
	Count =	[0] * len(DNA)

	for i in range (0, len(DNA)-k):
		PATTERN = Text(i, k, DNA)
		Count[i] = PatternCount(DNA, PATTERN)
			
	maxCount = max(Count)
	
	for i in range (0, len(DNA)-k):
		if (Count[i]) == maxCount:
			FrequentPatterns.append(Text(i, k, DNA))
	
	FrequentPatterns = remove_duplicates(FrequentPatterns)
		
	return FrequentPatterns

def ComputingFrequencies (DNA, k):
	FrequencyArray = [0] * (4**k)	
	
	for i in range (0, (len(DNA)-k)+1):
		pattern = Text(i, k, DNA)				
		j = PatternToNumber(pattern)
		FrequencyArray[j] = FrequencyArray[j] + 1
	
	return FrequencyArray

def ClumpFinding (GENOME, k, t, L):	
	FrequentPatterns = []
	clump = [0] * (4**k)
		
	for i in range (0, len(GENOME)-L+1):		
		text = GENOME[i:i+L]				
		FrequencyArray = ComputingFrequencies(text, k)
				
		for index in range (0, 4**k):
			if FrequencyArray[index] >= t:
				clump[index] = 1
	
	for index in range (0, 4**k):
		if clump[index] == 1:
			pattern = NumberToPattern(index, k)			
			FrequentPatterns.append(pattern)
		
	return FrequentPatterns
	
def BetterClumpFinding (GENOME, k, t, L):	
	FrequentPatterns = []
	clump = [0] * (4**k)

	text = GENOME[0:L]
	FrequencyArray = ComputingFrequencies(text, k)
	
	for index in range (0, 4**k):
		if FrequencyArray[index] >= t:
			clump[index] = 1
	
	for i in range (1, len(GENOME)-L+1):
		FirstPattern = GENOME[i-1:(i-1)+k]
		index = PatternToNumber(FirstPattern)
		FrequencyArray[index] = FrequencyArray[index] - 1
		
		LastPattern = GENOME[i+L-k:(i+L-k)+k]	
		index = PatternToNumber(LastPattern)
		FrequencyArray[index] = FrequencyArray[index] + 1
		
		if FrequencyArray[index] >= t:
			clump[index] = 1	
	
	for index in range (0, 4**k):
		if clump[index] == 1:
			pattern = NumberToPattern(index, k)			
			FrequentPatterns.append(pattern)
		
	return FrequentPatterns

"""
Following code is part of WEEK-2 of 
course BIO-INFORMATICS at COURSERA
"""
	
def Suffix(PATTERN):	
	return PATTERN[1:len(PATTERN)]

def Prefix(PATTERN):	
	return PATTERN[0:len(PATTERN)-1]
	
def FirstSymbol(PATTERN):
	return PATTERN[0]
	
def Skew(GENOME):

	skew = [0] * (len(GENOME)+1)
	
	minValue = skew[0]
		
	for i in range (0, len(GENOME)):		
		if GENOME[i] == "C":
			skew[i+1] = skew[i] - 1
		elif GENOME[i] == "G":
			skew[i+1] = skew[i] + 1
		elif GENOME[i] == "A" or GENOME[i] == "T":
			skew[i+1] = skew[i]
	
		if (skew[i+1] <= minValue):
			minValue = skew[i+1]
			
	# print ("min value: %d" %minValue)
	# print("Skew: ", skew) 
	"""E. Coli Skews: 3923620, 3923621, 3923622, 3923623"""
	
	# plt.plot(skew)
	# plt.ylabel('Skew')
	# plt.xlabel('Position')
	# plt.show()

	
	min_indexes = lambda x, xs: [i for (y, i) in zip(xs, range(len(xs))) if x == y]	
	indexes = min_indexes(minValue, skew)
	
	return indexes

def HammingDistance(p, q):
		
	dist = 0
	
	for i in range (0, len(p)):
		if p[i] != q[i]:
			dist = dist + 1
	
	return dist
	
def ApproximatePatternMatching(PATTERN, GENOME, d):
	
	indexes = []
	
	for i in range(0, len(GENOME)-len(PATTERN)+1):
		text = Text(i, len(PATTERN), GENOME)		
		
		if HammingDistance(PATTERN, text) <= d:
			indexes.append(i)
	
	return indexes
	
def ImmediateNeighbors (PATTERN):
	Neighborhood = []
	Neighborhood.append(PATTERN)
		
	for i in range (0, len(PATTERN)):
		symbol = PATTERN[i]
		
		SET = ['A', 'C', 'G', 'T']	
		SET.remove(symbol)
		
		for x in SET:
			Neighbor = PATTERN.replace(symbol, x)
			Neighborhood.append(Neighbor)
	
	return Neighborhood
		
def IterativeNeighbors (PATTERN, d):
	Neighborhood = []
	Neighborhood.append(PATTERN)
		
	for i in range(0, d):
		for string in Neighborhood:
			Neighborhood.extend(ImmediateNeighbors(string))			
			Neighborhood = remove_duplicates(Neighborhood)
			
			print ("i: %d, string: %s, N: " %(i, string), Neighborhood)
			input((">"))
	
	return Neighborhood
	
def Neighbors (PATTERN, d):
	if d == 0:
		return PATTERN
	
	SET = ['A', 'C', 'G', 'T']
	
	if len(PATTERN) == 1:
		return SET
		
	Neighborhood = []
	SuffixNeighbors = Neighbors(Suffix(PATTERN), d)
		
	for text in SuffixNeighbors:
		if HammingDistance(Suffix(PATTERN), text) < d:			
			for x in SET:
				Neighborhood.append(x+text)				
		else:
			Neighborhood.append(FirstSymbol(PATTERN)+text)
	
	return Neighborhood
			
def FrequentWordsWithMismatches (DNA, k, d):

	FrequencyArray = [0] * (4**k)	
	FrequentPatterns = []
	
	for i in range(0, len(DNA)-k+1):
		PATTERN = Text(i, k, DNA)
		Neighborhood = Neighbors(PATTERN, d)
				
		if type(Neighborhood) == str:
			j = PatternToNumber(Neighborhood)			
			FrequencyArray[j] = FrequencyArray[j] + 1
		else:
			for approximatePattern in Neighborhood:
				j = PatternToNumber(approximatePattern)				
				FrequencyArray[j] = FrequencyArray[j] + 1
		
	maxCount = max(FrequencyArray)
		
	for index in range (0, 4**k):
		if FrequencyArray[index] == maxCount:
			pattern = NumberToPattern(index, k)			
			FrequentPatterns.append(pattern)
	
	return FrequentPatterns

def FrequentWordsWithMismatchesAndReverseComplement (DNA, k, d):

	FrequencyArray1 = [0] * (4**k)	
	FrequencyArray2 = [0] * (4**k)	
	FrequentPatterns = []
	
	for i in range(0, len(DNA)-k+1):
		PATTERN = Text(i, k, DNA)
		Neighborhood1 = Neighbors(PATTERN, d)
		Neighborhood2 = Neighbors(ReverseComplement(PATTERN), d)
				
		if type(Neighborhood1) == str:
			j = PatternToNumber(Neighborhood1)			
			FrequencyArray1[j] = FrequencyArray1[j] + 1
		else:
			for approximatePattern in Neighborhood1:
				j = PatternToNumber(approximatePattern)				
				FrequencyArray1[j] = FrequencyArray1[j] + 1
		
		
		if type(Neighborhood2) == str:
			j = PatternToNumber(Neighborhood2)			
			FrequencyArray2[j] = FrequencyArray2[j] + 1
		else:
			for approximatePattern in Neighborhood2:
				j = PatternToNumber(approximatePattern)				
				FrequencyArray2[j] = FrequencyArray2[j] + 1
				
				
	maxCount = max(np.add(FrequencyArray1, FrequencyArray2))
	print ("\nMax Count: %d.\n" %maxCount)
	for index in range (0, 4**k):
		if FrequencyArray1[index]+FrequencyArray2[index] == maxCount:
			pattern = NumberToPattern(index, k)			
			FrequentPatterns.append(pattern)
	
	return FrequentPatterns	

def FrequentWordsWithMismatchesBySorting (DNA, k, d):

	FrequentPatterns = []
	Neighborhoods = []
					
	for i in range(0, len(DNA)-k+1):
		PATTERN = Text(i, k, DNA)
		Neighborhood = Neighbors(PATTERN, d)
		
		for approximatePattern in Neighborhood:
			j = PatternToNumber(approximatePattern)
			FrequencyArray[j] = FrequencyArray[j] + 1
		
	maxCount = max(FrequencyArray)
	
	for index in range (0, 4**k):
		if FrequencyArray[index] == maxCount:
			pattern = NumberToPattern(index, k)			
			FrequentPatterns.append(pattern)
	
	return FrequentPatterns

"""
Following code is part of WEEK-3 of 
course BIO-INFORMATICS at COURSERA
"""

def MotifEnumeration(DNA, k, d):
	PATTERNS = []
	
	for pat in DNA:
		# print("PAT:========> %s" %pat)
		for i in range(0, len(pat)-k+1):
			kmer = Text(i, k, pat)
			neighbors = Neighbors(kmer, d)
						
			# print("kmer: %s, N-L: %d,\n Neighbors: " %(kmer, len(neighbors)), neighbors)			
			# print(type(neighbors))
			# input((">"))
						
			if type(neighbors) == str:									
				flag = 1
				# print("n: %s" %neighbors)
				for pat2 in DNA:
					# print("pat2: %s" %pat2)
					aa = ApproximatePatternMatching(neighbors, pat2, d)					
					if len(aa) == 0:
						flag = 0
						break
				
				if flag == 1:
					PATTERNS.append(neighbors)
			else:				
				for n in neighbors:
					flag = 1
					# print("n: %s" %n)
					for pat2 in DNA:
						# print("pat2: %s" %pat2)
						aa = ApproximatePatternMatching(n, pat2, d)
						
						if len(aa) == 0:
							flag = 0
							break
							
					if flag == 1:
						PATTERNS.append(n)		
				
			# print("Patterns: ", PATTERNS)
			# print("\n\n")
			
	PATTERNS = remove_duplicates(PATTERNS)
	PATTERNS.sort()
	return PATTERNS
				
def DistanceBetweenPatternAndStrings(PATTERN, DNA):
	k = len(PATTERN)
	distance = 0
	
	for string in DNA:		
		hammingDistance = 100000000;
		for i in range(0, len(string)-len(PATTERN)+1):
			PATTERN_BAR = Text(i, k, string)
			HD = HammingDistance(PATTERN, PATTERN_BAR)
			
			if hammingDistance > HD:			
				hammingDistance = HD		
		
		distance = distance + hammingDistance
		
	return distance
	
def MedianString(DNA, k):
	distance = 100000000
	Median = 0
	
	for i in range (0, (4**k)):
		PATTERN = NumberToPattern(i, k)
		currentDistance = DistanceBetweenPatternAndStrings(PATTERN, DNA)
		if distance > currentDistance:
			distance = currentDistance
			Median = PATTERN
	
	return Median

def PatternProbabilityFromProfile (PATTERN, PROFILE):
	
	code = {'A': 0, 'C': 1, 'G': 2, 'T': 3}	
	prob = 1
	
	for c in range (0, len(PATTERN)):
		
		row = code[PATTERN[c]]						
		prob = prob*PROFILE[row][c]
		
		if prob == 0:
			return 0		
	
	return prob

def ProfileMostProbableKMER (string, k, PROFILE):
	
	probableKMER = []
	highestProbableKMER = 0
	highestProb = -10000000000
	
	for i in range (0, len(string)-k+1):
		kmer = Text(i, k, string)
		prob = PatternProbabilityFromProfile(kmer, PROFILE)
			
		if highestProb < prob:
			highestProb = prob
			highestProbableKMER = kmer
	
	return highestProbableKMER
		
def FindConsensus(PROFILE):
	
	code = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}	
	consensus = 0
	
	pr = np.array(PROFILE)
	max = pr.argmax(axis=0)
		
	consensus = code[max[0]]
	for i in range(1, len(max)):
		consensus = consensus + code[max[i]]
		
	return consensus
	
def Score(MOTIFS, PROFILE):
	
	consensus = FindConsensus(PROFILE)
	score = 0
	
	for m in MOTIFS:
		m = ''.join(m)
		dist = HammingDistance(consensus, m)
		score = score + dist
	
	
	return score

def CreateProfileMatrix(MOTIFS):
	
	t = len(MOTIFS)
	k = len(MOTIFS[0])	
	
	code = {'A': 0, 'C': 1, 'G': 2, 'T': 3}	
	MATRIX = [[0.0]*k for _ in range(4)]
		
	for c in range (0, k):		
		for r in range(0, t):				
			MATRIX[code[MOTIFS[r][c]]][c] = (MATRIX[code[MOTIFS[r][c]]][c] + 1)/t

	return MATRIX

def CreateLaplaceProfileMatrix(MOTIFS, laplace_const):
	
	t = len(MOTIFS)
	k = len(MOTIFS[0])	
	
	code = {'A': 0, 'C': 1, 'G': 2, 'T': 3}	
	
	MATRIX = [[1/(4+laplace_const)]*k for _ in range(4)]
	
	for col in range (0, k):		
		for row in range(0, t):				
			MATRIX[code[MOTIFS[row][col]]][col] = MATRIX[code[MOTIFS[row][col]]][col] + 1/(4+laplace_const)
	
	return MATRIX
		
def	GreedyMotifSearch(DNA, k, t):
	
	BestMotifs = []
	BestProfile = []
		
	for i in range (0, t):
		row = list(DNA[i])		
		BestMotifs.append(row[0:k])
		BestProfile = CreateProfileMatrix(BestMotifs)
			
	for i in range (0, len(DNA[0])-k+1):
		motifs_str = []
		motifs = []
		
		motifs.append(list(Text(i, k, DNA[0])))
		motifs_str.append(Text(i, k, DNA[0]))
						
		laplace_const = 1		
		for j in range(1, t):
			PROFILE = CreateLaplaceProfileMatrix(motifs, laplace_const)
			# PROFILE = CreateProfileMatrix(motifs)
			laplace_const += 1
			nextMotif = ProfileMostProbableKMER(DNA[j], k, PROFILE)
			motifs.append(list(nextMotif))
			motifs_str.append(nextMotif)
		
		if Score(motifs, PROFILE) < Score(BestMotifs, BestProfile):
			BestMotifs = motifs_str
			BestProfile = PROFILE		
	
	return BestMotifs

"""
Following code is part of WEEK-4 of 
course BIO-INFORMATICS at COURSERA
"""

def SelectRandomKMERS(DNA, k):
	rows = len(DNA)
	cols = len(DNA[0])
	
	kmers = []
	for i in range (0, rows):
		r = random.randint(0, cols-k)
		kmers.append(list(Text(r, k, DNA[i])))
			
	return kmers	

def FindMotifs(Profile, DNA, k):
	
	Motifs = []
	
	for i in range(0, len(DNA)):
		Motifs.append(list(ProfileMostProbableKMER(DNA[i], k, Profile)))
	
	return Motifs
	
def RandomizedMotifSearch(DNA, k, t):
		
	Motifs = SelectRandomKMERS(DNA, k)
	BestMotifs = Motifs
	BestProfile = CreateProfileMatrix(BestMotifs)
	BestScore = Score(BestMotifs, BestProfile)
	
	laplace_const = 1		
	while True:
		Profile = CreateLaplaceProfileMatrix(Motifs, laplace_const)
		Motifs = FindMotifs(Profile, DNA, k)
		laplace_const += 1
		sc = Score(Motifs, Profile)
			
		if sc < Score(BestMotifs, BestProfile):
			BestMotifs = Motifs
			BestProfile = Profile
			BestScore = sc
		else:
			return (BestMotifs, BestProfile, BestScore)			

def PatternProbabilityFromProfile2 (PATTERN, PROFILE):
	
	code = {'A': 0, 'C': 1, 'G': 2, 'T': 3}	
	prob = 1
	
	for c in range (0, len(PATTERN)):
		
		row = code[PATTERN[c]]						
		prob = prob*PROFILE[row][c]
		
		if prob == 0:
			return 0		
	
	return prob

def ProfileMostProbableKMER2 (string, k, PROFILE):
	
	probableKMER = []
	highestProbableKMER = 0
	highestProb = -10000000000
	
	sum = 0
	kmer = []
	prob = []
	
	for i in range (0, len(string)-k+1):
		kmer.append(Text(i, k, string))
		prob.append(PatternProbabilityFromProfile(kmer[i], PROFILE))
		sum += prob[i]
			
	pr = np.array(prob)
	pr = np.divide(pr, sum)
	prob = pr.tolist()
	
	for l in range (0, len(prob)):
		if highestProb < prob[l]:
			highestProb = prob[l]
			highestProbableKMER = kmer[l]
	
	return highestProbableKMER
			
def CreateLaplaceProfileMatrix2(MOTIFS):
	
	t = len(MOTIFS)
	k = len(MOTIFS[0])	
	
	code = {'A': 0, 'C': 1, 'G': 2, 'T': 3}	
	
	# print("==========", k, t, laplace_const)
	
	MATRIX = [[1]*k for _ in range(4)]
	
	# print(MATRIX)
	# print(MOTIFS)
	
	col_sum = 0
	for col in range (0, k):		
		for row in range(0, t):				
			MATRIX[code[MOTIFS[row][col]]][col] = MATRIX[code[MOTIFS[row][col]]][col] + 1					
	
	
	# for row in range(0, 4):				
		# col_sum += MATRIX[row][0]
	
	# print("\n")
	# print(MATRIX)
	# print("COL SUM: %d " %col_sum)
	
	# M = np.array(MATRIX)
	# M = np.divide(M, col_sum)
	# MATRIX = M.tolist()	
	
	# print(MATRIX)	
	# input((">"))
	
	return MATRIX
			
def GibbsSampler(DNA, k, t, N):
	
	BestMotifs = [] 
	BestProfile = []
		
	limit = 20
	minScore = 1000000000000
	
	for j in range (0, N):
					
		if j < limit:	
			random.seed()
			motifs_temp, profiles_temp, sc = RandomizedMotifSearch(DNA, k, t)
								
			if sc < minScore:
				minScore  = sc
				Motifs = motifs_temp
				Profile = profiles_temp		
				BestMotifs = Motifs
				BestProfile = Profile
			
			print("--", j, sc, minScore)	
		else:
			# print("After Randomization:")
			# print("Min Score: %d" %minScore)
			# input((">"))		
			
			i = random.randint(0, t-1)
			del Motifs[i]
			Profile = CreateLaplaceProfileMatrix2(Motifs)
			Motif_i = ProfileMostProbableKMER2(DNA[i], k, Profile)
			Motifs.insert(i, list(Motif_i))
		
			# print("Random I: %d" %i)
			# print(Profile)
			# print(Motif_i)	
			#
			
			sc = Score(Motifs, Profile)
			# bs = Score(BestMotifs, BestProfile)
					
			if sc < minScore:
				BestMotifs = Motifs
				BestProfile = Profile
				minScore = sc
	
			print(j, sc, minScore)	
	
	return BestMotifs

#______________________________________________________________________________		

"""
GENOME SEQUENCING: WEEK 01
"""

def Composition (DNA, k):
	
	kmers = []	
	for i in range (0, len(DNA[0])-k+1):
		kmers.append(Text(i, k, DNA[0]))
		
	kmers.sort()
	return kmers

def GenomePathString (kmers):
	genome = kmers[0]
		
	for j in range (1, len(kmers)):		
		genome += (kmers[j])[-1]
	
	return genome
		
def OverlapGraph (kmers):
	connectivityGraph = []
	
	prefix = []
	suffix = []
	
	for i in range(0, len(kmers)):
		suffix.append(Suffix(kmers[i]))
		prefix.append(Prefix(kmers[i]))
		
	# print(suffix)
	# print(prefix)		
	
	indices = []
	for i in range(0, len(kmers)):
		if suffix[i] in prefix:
			indices.append(prefix.index(suffix[i]))
		else:
			indices.append(-1)
	
	return indices

def DeBruijnGraph (k, string):
	
	connectionDict = {}
	
	for i in range (0, len(string) - k + 2):
		
		kmer_now = Text(i, k-1, string)
		kmer_next = Text(i+1, k-1, string)
		
		# print(kmer_now, kmer_next)
		
		if kmer_now not in connectionDict and len(kmer_next) == k-1:
			connectionDict[kmer_now] = kmer_next
		elif len(kmer_next) == k-1:
			connectionDict[kmer_now] = connectionDict[kmer_now] + "," + kmer_next
			temp = connectionDict[kmer_now].split(',')			
			connectionDict[kmer_now] = ','.join(sorted(temp))
		
	connectionDict = collections.OrderedDict(sorted(connectionDict.items()))
	# print("\n")
	# print(connectionDict)
	
	return connectionDict

def DeBruijnGraphFromKmers (kmers):
	
	connectionDict = {}
	connectionArray = {}
			
	for i in range (0, len(kmers)):
		prefix = Prefix(kmers[i])
		suffix = Suffix(kmers[i])
		
		connectionArray[prefix] = prefix		
		connectionArray[suffix] = suffix

	for i in range(0, len(kmers)):
		
		prefix = Prefix(kmers[i])
		suffix = Suffix(kmers[i])
		
		if prefix not in connectionDict:
			connectionDict[connectionArray[prefix]] = suffix
		else:
			connectionDict[connectionArray[prefix]] = connectionDict[connectionArray[prefix]] + "," + suffix
			temp = connectionDict[connectionArray[prefix]].split(',')			
			connectionDict[connectionArray[prefix]] = ','.join(sorted(temp))
		
	
	# connectionDict = collections.OrderedDict(sorted(connectionDict.items()))
	# print(connectionDict)		
	# print("\n")
	
	return connectionDict

"""
GENOME SEQUENCING: WEEK 02
"""

def getAdjacencyMatrix (deBruijnGraph):
	
	temp = []
	for key, value in deBruijnGraph.items():
		temp.append("%s->%s" %(key, value))
	
	
	ADJACENCY = collections.defaultdict(list)
	inEdges = collections.defaultdict(list)

	for line in temp:
		conn = re.findall('\w+|\d+', line)		
		# print(line, conn)
		for i in range (1, len(conn)):
			ADJACENCY[conn[0]].append(conn[i])
			inEdges[conn[i]].append(conn[0])
	
	return ADJACENCY, inEdges

def EulerianCycle (ADJACENCY, highestKey):

	STACK = []
	currentNode = highestKey
	CIRCUIT = []
	
	# print("\nSTACK: ", STACK)
	# print("Current Node: ", currentNode)
	# print("Circuit: ", CIRCUIT)
		
	while True:		
		if currentNode in ADJACENCY and len(ADJACENCY[currentNode]) != 0:
			# Add the currentNode to the stack, take any of its neighbors, 
			# remove the edge between selected neighbor and that vertex, 
			# and set that neighbor as the current vertex.
						
			STACK.append(currentNode)
			neighborList = ADJACENCY[currentNode]
			neighbor = neighborList.pop()	
			
			# print("\nCN: ", currentNode)		
			# print(neighborList, neighbor, len(ADJACENCY[currentNode]), len(neighborList))
			# input((">"))
						
			#Delete edge between currentNode and the selected neighbor
			for index, key in enumerate(ADJACENCY[currentNode]):
				if key == neighbor:
					ADJACENCY[currentNode].pop(index)
					
			currentNode = neighbor			
		else:					
			# If currentNode has no neighbors, then add it add it to circuit, 
			# remove the last Node from the stack and set it as the currentNode
			CIRCUIT.append(str(currentNode))
			
			if len(STACK) == 0:
				break
			else:
				currentNode = STACK.pop()	
		
		# print("\nSTACK: ", STACK)
		# print("Current Node: ", currentNode)
		# print(ADJACENCY)
		# print("Circuit: ", CIRCUIT)
		# input((">"))	
		
	CIRCUIT = CIRCUIT[::-1] 		
	
	return CIRCUIT

def EulerianPath (ADJACENCY, startNode):
			
	STACK = []
	currentNode = startNode
	CIRCUIT = []
			
	while True:	
		# print("\nCN: ", currentNode)
		# input((">"))
		if currentNode in ADJACENCY and len(ADJACENCY[currentNode]) != 0:
			# Add the currentNode to the stack, take any of its neighbors, 
			# remove the edge between selected neighbor and that vertex, 
			# and set that neighbor as the current vertex.
						
			STACK.append(currentNode)
			neighborList = ADJACENCY[currentNode]					
			neighbor = neighborList.pop()			
			
			# print(neighborList, neighbor, len(ADJACENCY[currentNode]), len(neighborList))			
									
			#Delete edge between currentNode and the selected neighbor
			for index, key in enumerate(ADJACENCY[currentNode]):
				if key == neighbor:
					ADJACENCY[currentNode].pop(index)
				
			currentNode = neighbor			
		else:					
			# If currentNode has no neighbors, then add it add it to circuit, 
			# remove the last Node from the stack and set it as the currentNode
			CIRCUIT.append(str(currentNode))
			# print("STACK: ", STACK)
			# print("CIRCUIT: ", CIRCUIT)
			
			if len(STACK) == 0:
				break
			else:
				currentNode = STACK.pop()	
		
		
		# print("ADJACENCY: ",ADJACENCY)
		
		
	CIRCUIT = CIRCUIT[::-1] 	
	
	return CIRCUIT
	
def kUniversalString (k):
	string = []
	
	kmers = [ ''.join(x) for x in product('01', repeat=k) ]
	
	deBruijn = DeBruijnGraphFromKmers(kmers)
	ADJACENCY = getAdjacencyMatrix(deBruijn)
	
	# print(deBruijn)
	# print(ADJACENCY)
	
	highestKey = -100000
	highestValue = -100000
	for items in sorted(ADJACENCY.keys()):  
		# print("%s: Fan-out: %d" %(items, len(ADJACENCY[items]))) 			
		if len(ADJACENCY[items]) > highestValue:
			highestKey = items
			highestValue = len(ADJACENCY[items])
	
	
	eulerCycle = EulerianCycle(ADJACENCY, highestKey)
	eulerCycle[0] = kmers[0]		
	string = GenomePathString(eulerCycle)			
	string = string[0:-k]
	
	return string

def StringSpelledByGappedPatterns (GappedPatterns, k, d):

	FirstPatterns = []
	SecondPatterns = []
		
	for i in range (0, len(GappedPatterns)):
		pattern = GappedPatterns[i]
		pattern = pattern.split('|')
		FirstPatterns.append(pattern[0])
		SecondPatterns.append(pattern[1])
			
	# print(FirstPatterns)
	# print(SecondPatterns)
	
	PrefixString = GenomePathString (FirstPatterns)
	SuffixString = GenomePathString (SecondPatterns)
	
	# print(PrefixString)
	# print(SuffixString)
		
	for i in range (k+d, len(PrefixString)):
		if PrefixString[i] != SuffixString[i-k-d]:
			return "there is no string spelled by the gapped patterns"
		
	string = PrefixString + SuffixString[-(k+d):]
	return string

def DeBruijnGraphFromPairedKmers (kmers):
	
	connectionDict = {}
	connectionArray = {}
	
	# print(kmers)
			
	for i in range (0, len(kmers)):
		patterns = kmers[i].split('|')
		
		prefix_FP = Prefix(patterns[0])
		suffix_FP = Suffix(patterns[0])
		
		prefix_SP = Prefix(patterns[1])
		suffix_SP = Suffix(patterns[1])
		
		prefix = prefix_FP + "|" + prefix_SP
		suffix = suffix_FP + "|" + suffix_SP
		
		connectionArray[prefix] = prefix		
		connectionArray[suffix] = suffix
	
	# print(connectionArray)
	# input((">"))

	for i in range(0, len(kmers)):
		
		patterns = kmers[i].split('|')
				
		prefix_FP = Prefix(patterns[0])
		suffix_FP = Suffix(patterns[0])
		
		prefix_SP = Prefix(patterns[1])
		suffix_SP = Suffix(patterns[1])
		
		prefix = prefix_FP + "|" + prefix_SP
		suffix = suffix_FP + "|" + suffix_SP
		
		# print(prefix, suffix)
		
		if prefix not in connectionDict:
			connectionDict[connectionArray[prefix]] = suffix
		else:			
			connectionDict[connectionArray[prefix]] = connectionDict[connectionArray[prefix]] + "," + suffix
			temp = connectionDict[connectionArray[prefix]].split(',')			
			connectionDict[connectionArray[prefix]] = ','.join(sorted(temp))
		
	
	# connectionDict = collections.OrderedDict(sorted(connectionDict.items()))
	# print(connectionDict)		
	# print("\n")
	
	return connectionDict
	
def getAdjacencyMatrixForPairedReads (pairedDeBruijnGraph):
	
	temp = []
	for key, value in pairedDeBruijnGraph.items():
		temp.append("%s->%s" %(key, value))
	
	
	ADJACENCY = collections.defaultdict(list)
	inEdges = collections.defaultdict(list)

	for line in temp:
		conn = re.findall('\w+|\d+', line)		
		# print(line, conn)
		
		for i in range (2, len(conn), 2):			
			str = conn[0] + "|" + conn[1]
			str2 = conn[i] + "|" + conn[i+1]
		
			ADJACENCY[str].append(str2)
			inEdges[str2].append(str)
			
	return ADJACENCY, inEdges
	
def StringReconstructionFromReadPairs (kmers, k, d):

	pairedDeBruijnGraph = DeBruijnGraphFromPairedKmers(kmers)	
	ADJACENCY, inEdges  = getAdjacencyMatrixForPairedReads(pairedDeBruijnGraph)
	
	temp = {**ADJACENCY, **inEdges}	
	for items in sorted(temp.keys()):						
		if items in inEdges and items not in ADJACENCY:
			len(ADJACENCY[items]) # This inserts a dummy edge with zero output
	
		fanIns = len(inEdges[items])
		fanOuts = len(ADJACENCY[items])
		
		if fanOuts - fanIns == 1:
			startNode = items
		elif fanIns - fanOuts == 1:
			endNode = items
	
	eulerPath = EulerianPath(ADJACENCY, startNode)			
	genomeString = StringSpelledByGappedPatterns(eulerPath, k, d)


	return genomeString

def MaximalNonBranchingPaths (Graph):

	fanIns	=	collections.defaultdict(list)
	fanOuts	=	collections.defaultdict(list)

	nodes = []
	
	for i in range (len(Graph)):
		line = re.findall('\w+|\d+', Graph[i])	
		
		if line[0] not in nodes:
			nodes.append(line[0])
		
		fanIns[line[0]]
		
		for k in range (1, len(line)):
			fanOuts[line[0]].append(line[k])
			fanIns[line[k]].append(line[0])
			
			fanOuts[line[k]]
			fanIns[line[k]]
			
			if line[k] not in nodes: 
				nodes.append(line[k])	
	
		
	Paths = []
	NonBranchingPath = []	
	visited = {}
	
	for i in range (len(nodes)+1):
		visited[str(i)] = False
		
	#Explore each component alone		
	for node in (nodes):
		# print("+++",node, len(fanIns[node]), len(fanOuts[node]))		
		
		if len(fanIns[node]) != 1 or len(fanOuts[node]) != 1:			
			if len(fanOuts[node]) > 0:
				
				visited[node] = True				
				for edge in (fanOuts[node]):
					# print(node, edge)
					NonBranchingPath.append(node)
					NonBranchingPath.append(edge)
					visited[edge] = True
					
					while(len(fanIns[edge]) == 1 and len(fanOuts[edge]) == 1):											
						nextNode = fanOuts[edge].pop()
						# print("--", edge, nextNode)	
						edge = nextNode	
						visited[edge] = True						
						NonBranchingPath.append(edge)
						
					
					Paths.append(NonBranchingPath)
					NonBranchingPath = []	
		
	
	nodes = []
	for i in (visited):
		if visited[i] == False:			
			nodes.append(i)
	
	NonBranchingPath = []
	
	for node in (nodes):		
		# if len(fanIns[node]) != 1 or len(fanOuts[node]) != 1:			
		if len(fanOuts[node]) > 0:			
			
			visited[node] = True				
			
			for edge in (fanOuts[node]):
				# print(node, edge)
				NonBranchingPath.append(node)
				NonBranchingPath.append(edge)
				visited[edge] = True
				
				while(len(fanIns[edge]) == 1 and len(fanOuts[edge]) == 1):											
					nextNode = fanOuts[edge].pop()
					# print("--", edge, nextNode)	
					edge = nextNode

					if visited[edge] == True:
						NonBranchingPath.append(edge)
						break
						
					visited[edge] = True						
					NonBranchingPath.append(edge)
					
				
				Paths.append(NonBranchingPath)
				NonBranchingPath = []	
		
	# print("FanOuts: ", fanOuts)
	# print("FanIns: ", fanIns)
	# print("Nodes: ", nodes)
			
	return Paths

def ContigsFromBranchingPaths (Paths):
	
	contigs = []
	
	for text in (Paths):
		temp = text[0]
		for k in range(1, len(text)):
			temp = temp + (text[k])[-1]
		
		contigs.append(temp)
		
		
	return sorted(contigs)

"""
GENOME SEQUENCING: WEEK 03
"""

def NonOverLappingKmers (Genome, k):
	
	kmers = []			
	for i in range (0, len(Genome), k):
		kmers.append(Text(i, k, Genome))
		
	return kmers	

def OverLappingKmers (Genome, k):
	
	kmers = []			
	for i in range (0, len(Genome)-k+1):
		kmers.append(Text(i, k, Genome))
		
	return kmers	
	
def RNA2AminoString (RNA):

	txt = open('rna-codon.txt', 'r')
	
	codons = collections.defaultdict()	
	for line in txt:	
		temp = re.findall('\w+|\d+|\*', line)				
		codons[temp[0]] = temp[1]
	txt.close()
	
	kmers = NonOverLappingKmers(RNA, 3)
	
	peptide = []
	for i in (kmers):
		peptide.append(codons[i])
	
	peptide = ''.join(peptide)
	
	return peptide

def DNA_TO_RNA (DNA):
	
	RNA = []	
	for i in range(len(DNA)):
		if DNA[i] == "T":
			RNA.append("U")
		else:
			RNA.append(DNA[i])			
			
	RNA = ''.join(RNA)
	
	return RNA

def RNA_TO_DNA (RNA):
	
	DNA = []	
	for i in range(len(RNA)):
		if RNA[i] == "U":
			DNA.append("T")
		else:
			DNA.append(RNA[i])			
			
	DNA = ''.join(DNA)	
	return DNA
	
def Peptide_Encoding (DNA, Peptide):

	k = len(Peptide)*3
	
	txt = open('rna-codon.txt', 'r')	
	codons = collections.defaultdict()	
	for line in txt:	
		temp = re.findall('\w+|\d+|\*', line)				
		codons[temp[0]] = temp[1]
	txt.close()
	
	kmers = OverLappingKmers(DNA, k)
			
	peptide_encoders = []	
	for kmer in (kmers):
		RC = ReverseComplement(kmer)		
		RNA_kmer = DNA_TO_RNA(kmer)
		RC_kmer = DNA_TO_RNA(RC)
		
		temp1 = ""
		temp2 = ""
		
		for m in range (0, len(kmer), 3):			
			temp1 += codons[RNA_kmer[m:m+3]] 			
			temp2 += codons[RC_kmer[m:m+3]] 				
		
		if temp1 == Peptide or temp2 == Peptide:
			peptide_encoders.append(kmer)

	
	return peptide_encoders

def LinearSpectrum (Peptide, AminoAcid):

	PrefixMass = []
	PrefixMass.append(0)
	
	for i in range (1, len(Peptide)+1):
		PrefixMass.append(PrefixMass[i-1] + AminoAcid[Peptide[i-1]])
	
	Lin_Spectrum = []
	Lin_Spectrum.append(0)	
	for i in range (len(Peptide)):
		for j in range (i+1, len(Peptide)+1):
			temp = PrefixMass[j] - PrefixMass[i]			
			Lin_Spectrum.append(temp)
    
	return sorted(Lin_Spectrum)	
	
def Cyclospectrum (Peptide, AminoAcid):

	PrefixMass = []
	PrefixMass.append(0)
	
	for i in range (1, len(Peptide)+1):
		PrefixMass.append(PrefixMass[i-1] + AminoAcid[Peptide[i-1]])
	
	
	peptideMass = PrefixMass[-1]
	Cyclic_Spectrum = []
	Cyclic_Spectrum.append(0)	
	
	# print(PrefixMass, peptideMass)
	
	for i in range (len(Peptide)):
		for j in range (i+1, len(Peptide)+1):	
			temp = PrefixMass[j] - PrefixMass[i]
			# print(temp)
			Cyclic_Spectrum.append(temp)
			
			if i > 0 and j < len(Peptide):
				temp = peptideMass - (PrefixMass[j] - PrefixMass[i])
				# print("--", temp)
				Cyclic_Spectrum.append(temp)
    
	return sorted(Cyclic_Spectrum)
	
def Expand (Peptides, AminoAcid):

	Peptides_temp = []		
	
	if len(Peptides) == 1 and Peptides[0] == "X":
		Peptides = []
		for key in AminoAcid.keys():			
			Peptides.append(key)
		return Peptides
	else:		
		for pep in (Peptides):
			for key in AminoAcid.keys():				
				Peptides_temp.append(pep+key)			
		
	return Peptides_temp

def ComputeMass (Peptide, AminoAcid):
	mass = 0
	pep = list(Peptide)
	
	for i in (pep):
		mass += AminoAcid[i]
		
	return mass

def CompareLinearSpectrumWithSpectrum (LinSpectrum, Spectrum):
	
	# print(LinSpectrum, Spectrum)
	# input((">"))
	
	for i in (LinSpectrum):
		# print(i, LinSpectrum.count(i), Spectrum.count(i))
		if LinSpectrum.count(i) > Spectrum.count(i):			
			# input((">>>>"))
			return False
	
	
	return True
		
def CyclopeptideSequencing (Spectrum, AminoAcid):
	
	Peptides = []
	result = []
	Peptides.append("X")
	
	Spectrum = sorted(Spectrum)
	parentMass = Spectrum[-1]
				
	while len(Peptides) > 0:
		Peptides = Expand(Peptides, AminoAcid)		
		Peptides_temp = []
		for pep in (Peptides):
			mass = ComputeMass(pep, AminoAcid)			
						
			if mass in Spectrum:
				Peptides_temp.append(pep)
				if mass == parentMass:
					LinSpectrum = LinearSpectrum(pep, AminoAcid)
					if (CompareLinearSpectrumWithSpectrum(LinSpectrum, Spectrum)):												
						result.append(pep)							
			
		Peptides = list(Peptides_temp)		
	
	
	result = result[::-1]
	final = []
	for i in (result):
		temp = []
		temp.append(str(AminoAcid[i[0]]))		
		for k in range(1, len(i)):
			temp.append(str(AminoAcid[i[k]]))
		
		final.append('-'.join(temp))
	
	return (final)

	
"""
GENOME SEQUENCING: WEEK 04
"""
	
def CycloScoring (Peptide, Spectrum, AminoAcid):	
	
	if len(Peptide) == 1 and Peptide == " ":
		return 0
		
	cycloSpec = Cyclospectrum(Peptide, AminoAcid)
	spec = list(Spectrum)	
	
	# spec = [242, 242]
	# cycloSpec   = [242,  242, 242,  242]
	
	score = 0
	for mass in (cycloSpec):		
		if mass in spec:
			score += 1
			spec.remove(mass)
		
	return score

def LinearScoring (Peptide, Spectrum, AminoAcid):	
	
	if len(Peptide) == 1 and Peptide == " ":
		return 0
		
	linearSpec = LinearSpectrum(Peptide, AminoAcid)
	spec = list(Spectrum)	
	
	score = 0
	for mass in (linearSpec):
		if mass in spec:
			score += 1
			spec.remove(mass)
		
	return score

def Trim (Leaderboard, Spectrum, N, AminoAcid):
	
	LinearScores = collections.defaultdict()
	
	for peptide in (Leaderboard):
		temp = LinearScoring(peptide, Spectrum, AminoAcid)
		LinearScores[peptide] = temp		
		
	sorted_x = dict(sorted(LinearScores.items(), key=operator.itemgetter(1), reverse=True))
	
	k = 1
	nth_key = 0
	LB = []
	for key in sorted_x:		
		if (k > N):
			if sorted_x[key] == sorted_x[nth_key]:
				LB.append(key)				
			elif sorted_x[key] < sorted_x[nth_key]:
				return LB
		elif k == N:
			nth_key = key
			LB.append(key)
		else:
			LB.append(key)
			
		k += 1	
	
	return LB

def LeaderboardCyclopeptideSequencing (Spectrum, N, AminoAcid):

	Leaderboard = []
	Leaderboard.append("X")
	LeaderPeptide = " "
	LeaderPeptide_Score = 0
	
	Spectrum = sorted(Spectrum)
	parentMass = Spectrum[-1]
	print("Parent Mass: ", parentMass)
	
	LeadingLinearPeptides = []
			
	while len(Leaderboard) > 0:
		Leaderboard = Expand(Leaderboard, AminoAcid)
		
		Leaderboard_temp = list(Leaderboard)
		for Peptide in (Leaderboard_temp):
			mass = ComputeMass(Peptide, AminoAcid)
		
			if mass in Spectrum:
				if mass == parentMass:							
					sc1 = CycloScoring(Peptide, Spectrum, AminoAcid)
					sc2 = CycloScoring(LeaderPeptide, Spectrum, AminoAcid)
					
					print(Peptide, sc1)
					# if sc1 == 83:
					LeadingLinearPeptides.append(Peptide)
					
					if sc1 > sc2:					
						LeaderPeptide = Peptide
						LeaderPeptide_Score = sc1
						# print("+++", LeaderPeptide.encode('ascii', 'ignore'), len(LeaderPeptide), LeaderPeptide_Score)					
						
			elif mass not in Spectrum or mass > parentMass:
				Leaderboard.remove(Peptide)						
			
		# print("---", LeaderPeptide.encode('utf-8'), len(LeaderPeptide), LeaderPeptide_Score)
		print("---", LeaderPeptide, len(LeaderPeptide), LeaderPeptide_Score)	
						
		Leaderboard = Trim(Leaderboard, Spectrum, N, AminoAcid)
	
	
	# for i in (LeaderPeptide):
		# print(i.encode('ascii', 'ignore'))
	
	temp = []
	temp.append(str(AminoAcid[LeaderPeptide[0]]))		
	for k in range(1, len(LeaderPeptide)):
		temp.append(str(AminoAcid[LeaderPeptide[k]]))		
	
	LeaderPeptide_Mass = '-'.join(temp)
	print(LeaderPeptide)
	
	
	print(LeadingLinearPeptides, len(LeadingLinearPeptides))	
	LP_83 = []
	for peptide in (LeadingLinearPeptides):
		sc1 = CycloScoring(peptide, Spectrum, AminoAcid)
		
		temp = []
		temp.append(str(AminoAcid[peptide[0]]))		
		for k in range(1, len(peptide)):
			temp.append(str(AminoAcid[peptide[k]]))	
				
		LP_83.append('-'.join(temp))
		# print(peptide.encode('ascii', 'ignore'), sc1)
		print('-'.join(temp), sc1)	
		
	return LeaderPeptide_Mass
    
def SpectralConvolution (Spectrum):

	Spectrum_X = Spectrum.copy()
	Spectrum_Y = Spectrum.copy()
	
	spec_conv = []
	
	for i in Spectrum_X:
		for j in Spectrum_Y:
			temp = i - j
			if temp > 0:
				spec_conv.append(temp)
	
	return spec_conv
	
def ConvolutionCyclopeptideSequencing (M, N, conv_spectrum, Spectrum):

	freq = Counter(conv_spectrum)	
	sorted_keys = sorted(freq,  key=freq.__getitem__, reverse=True)	
	
	# print(freq)
	# print(sorted_keys)	

	temp = []
	for i in (sorted_keys):
		if i >= 57 and i <= 200:
			temp.append(i)			
			
	# print("TEMP: ", temp)
	M_Freq = []
	for i in range(0, len(temp)):
		if len(M_Freq) > M and freq[temp[i]] != freq[temp[i+1]]:
			M_Freq.append(temp[i])
			break
		else:
			M_Freq.append(temp[i])
	
	# print(M_Freq, len(M_Freq))
	
	
	AminoAcid = collections.defaultdict()
	for k in (M_Freq):
		AminoAcid[chr(k)] = k
		# print(k, chr(k))
		
	print(AminoAcid)
	
	LeaderPeptide = LeaderboardCyclopeptideSequencing(Spectrum, N, AminoAcid)
	
	return LeaderPeptide

	
###############################################################################
script, inputFile = argv
txt = open(inputFile, 'r')
###############################################################################

# DNA	=	txt.readline().rstrip('\n').split(" ")
# LB	=	txt.readline().rstrip('\n').split(" ")


M = txt.readline().rstrip("\n")
M = int(M)
N = txt.readline().rstrip("\n")
N = int(N)

k = txt.readline().rstrip("\n").split(" ")
# k, t, N = k.split(" ")
# k = int(k)
# t = int(t)

Spectrum = []
for i in (k):
	Spectrum.append(int(i))
Spectrum = sorted(Spectrum)


AminoAcid = collections.defaultdict() 
for i in range(57,201):
	AminoAcid[chr(i)] = int(i)

	
# txt2 = open('integer_mass_table.txt', 'r')
# for line in txt2:   
    # temp = re.findall('\w+|\d+|\*', line)
    # AminoAcid[temp[0]] = int(temp[1])

conv_spectrum = SpectralConvolution(Spectrum)
print(ConvolutionCyclopeptideSequencing(M, N, conv_spectrum, Spectrum))
# print(AminoAcid)
# print(LeaderboardCyclopeptideSequencing(Spectrum, N, AminoAcid))



# print(Trim(LB, Spectrum, N, AminoAcid))
# print(CycloScoring(N, Spectrum, AminoAcid))
	
# print ("-"*50)
# print ("DNA: %s, PATTERN: %s" %(DNA, PATTERN))
# print ("DNA Length: %d, PATTERN Length: %d" %(len(DNA), len(PATTERN)))
# print("PAT: ", PATTERN)
# print("DNA: ", DNA)
# print ("Motifs: ", MOTIFS)
# print ("DNA1 Length: %d, DNA2 Length: %d" %(len(DNA), len(DNA2)))
# print ("Index: %d" %index)
# print ("k: %d, t: %d, N: %d" %(k, t, N))
# print ("k: %d, d: %d" %(k, d))
# print ("k: %d, L: %d, t:%d" %(k, L, t))
# print ("-"*50)
# input((">"))


# print(CyclicSpectrum(PATTERN, AminoAcid))

# count = 0
# ans = []
# t0 = time.time()
# for i in (DNA):
	# # print(i)
	# result = Peptide_Encoding(i, Peptide)
	
	# if len(result) > 0:
		# count += 1
		# ans.append(result)	

# t1 = time.time()
# total = t1-t0

# print(ans, count, total)