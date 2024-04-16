
import time, os, psutil, sys
#import matplotlib.pyplot as plt
mismatchCost = dict()
mismatchCost["A"] = dict()
mismatchCost["C"] = dict()
mismatchCost["T"] = dict()
mismatchCost["G"] = dict()
mismatchCost["A"]["A"] = 0
mismatchCost["A"]["C"] = 110
mismatchCost["A"]["G"] = 48 
mismatchCost["A"]["T"] = 94

mismatchCost["C"]["A"] = 110
mismatchCost["C"]["C"] = 0
mismatchCost["C"]["G"] = 118 
mismatchCost["C"]["T"] = 48

mismatchCost["G"]["A"] = 48
mismatchCost["G"]["C"] = 118
mismatchCost["G"]["G"] = 0
mismatchCost["G"]["T"] = 110

mismatchCost["T"]["A"] = 94
mismatchCost["T"]["C"] = 48
mismatchCost["T"]["G"] = 110 
mismatchCost["T"]["T"] = 0
gapCost = 30
outputFile = ''

def stringGenerator(params, indices):
	res = params
	prev = res
	for i in indices:
		res = res[0:i+1] + prev + res[i+1:]
		prev = res
	return res

def outputFuntion(A, X, Y, a, b):
	String1 = ""
	String2 = ""
	i, j = a, b
	while i > 0 and j > 0: 
		if A[i][j] == A[i-1][j-1] + mismatchCost[X[i-1]][Y[j-1]]:
			String1 = String1 + X[i-1]
			String2 = String2 + Y[j-1]
			i = i-1
			j = j-1
		elif A[i][j] == gapCost + A[i-1][j]:
			String1 = String1 + X[i-1]
			String2 = String2 + "_"
			i = i-1
		elif A[i][j] == gapCost + A[i][j-1]:
			String1 = String1 + "_"
			String2 = String2 + Y[j-1]
			j = j-1
	while i > 0:
		String1 = String1 + X[i-1]
		String2 = String2 + "_"
		i = i-1
	while j > 0:
		String1 = String1 + "_"
		String2 = String2 + Y[j-1]
		j = j-1

	p = String1[::-1]
	q = String2[::-1]


	global outputFile
	outputFile = outputFile + str(float(A[a][b])) + '\n' + p + '\n' + q + '\n' 
	

def sequenceForMinimalAlignment(X,Y):
	m = len(X)
	n = len(Y)
	A = [[0 for c in range(n+1)] for r in range(m+1)]
	for i in range(m+1):
		A[i][0] = gapCost * i
	for j in range(n+1):
		A[0][j] = gapCost * j
	for j in range(1, n+1):
		for i in range(1, m+1):
			A[i][j] = min(mismatchCost[X[i-1]][Y[j-1]]+A[i-1][j-1], gapCost+A[i-1][j], gapCost+A[i][j-1])
	outputFuntion(A, X, Y, m, n)


if __name__ == "__main__":
	stime=time.time()
	with open(sys.argv[1]) as f:
		data = f.readlines()
	f.close()
	n = len(data)
	temp1 = data[0].rstrip('\r\n')
	lentemp1 = len(temp1)
	indices_s1 = []
	i = 1
	
	while len(data[i])<5:
		indices_s1.append(int(data[i]))
		i = i+1
	j = i-1

	temp2 = data[i].rstrip('\r\n')
	lentemp2 = len(temp2)
	i = i+1
	indices_s2 = []
	while i<n:
		indices_s2.append(int(data[i].rstrip('\r\n')))
		i = i+1
	
	s1 = stringGenerator(temp1,indices_s1)
	s2 = stringGenerator(temp2,indices_s2)
	
	sequenceForMinimalAlignment(s1,s2)

	
	fileOutput = open("output.txt", "w")
	outputFile = outputFile + str(format(((time.time() - stime) * 1000),".4f")) + '\n' + str(float(psutil.Process(os.getpid()).memory_info().rss /1024))+'\n'
	fileOutput.write(outputFile)
	fileOutput.close()




