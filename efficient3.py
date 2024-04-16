import time, os, psutil, sys
gapCost = 30
mismatchCost = {'AA': 0, 'AC':110, 'CA': 110, 'AG': 48, 'GA': 48, 'CC':0, 'AT': 94, 'TA': 94, 'CG':118, 'GC': 118, 'GG': 0, 'TG':110, 'GT':110, 'TT':0, 'CT': 48, 'TC':48}
    
def stringGenerator(params, indices):
	res = params
	prev = res
	for i in indices:
		res = res[0:i+1]+prev+res[i+1:]
		prev = res
	return res

def outputFuntion(A, X, Y, a, b):
	String1 = ""
	String2 = ""
	i, j = a, b
	while j > 0 and i > 0: 
		if A[i][j] == A[i-1][j-1] + mismatchCost[X[i-1]+Y[j-1]]:
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
	return p, q

def baseAlignment(X, Y):
	m = len(X)
	n = len(Y)
	A = [[0 for c in range(n+1)] for r in range(m+1)]
	for i in range(m+1):
		A[i][0] = gapCost * i
	for j in range(n+1):
		A[0][j] = gapCost * j
	for j in range(1, n+1):
		for i in range(1, m+1):
			A[i][j] = min(mismatchCost[X[i-1]+Y[j-1]]+A[i-1][j-1], gapCost+A[i-1][j], gapCost+A[i][j-1])
	p, q = outputFuntion(A, X, Y, m, n)	
	return p, q, A[m][n]


def sequenceForMinimalAlignment(X, Y):
	n = len(X)
	m = len(Y)
	B = [[0 for i in range(m+1)] for j in range(2)]
	for i in range(m+1):
		B[0][i] = i * gapCost
	for i in range(1, n+1):
		B[1][0] = B[0][0] + gapCost
		for j in range(1, m+1):

			B[1][j] = min(B[0][j-1] + mismatchCost[X[i-1] + Y[j-1]],
                            B[0][j] + gapCost,
                            B[1][j-1] + gapCost )
		for i in range(0, m+1):
			B[0][i] = B[1][i]

	return B[1]


def sequenceForBackwardSpaceAlignment(X, Y):
	n = len(X)
	m = len(Y)
	B = [[0 for i in range(m+1)] for j in range(2)]
	for i in range(m+1):
		B[0][i] = i * gapCost
	for i in range(1, n+1):
		B[1][0] = B[0][0] + gapCost
		for j in range(1, m+1):
			B[1][j] = min(mismatchCost[X[n-i]+Y[m-j]]+B[0][j-1] , gapCost + B[0][j], gapCost + B[1][j-1])
		for i in range(0, m+1):
			B[0][i] = B[1][i]

	return B[1]


def DAC(X, Y):
	n = len(X)
	m = len(Y)
	if m <= 2 or n <= 2:
		return baseAlignment(X, Y)
	else:
		leftPart = sequenceForMinimalAlignment(X[:n//2], Y)
		rightPart = sequenceForBackwardSpaceAlignment(X[n//2:], Y)

		median = [leftPart[j]+rightPart[m-j] for j in range(m+1)]
		minimalCut = median.index(min(median))

		leftPart,rightPart,median = [],[],[]

		leftResult = DAC(X[:n//2], Y[:minimalCut])
		rightResult = DAC(X[n//2:], Y[minimalCut:])

		return [leftResult[i] + rightResult[i] for i in range(3)]


if __name__ == "__main__":
	stime = time.time()
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

	r = DAC(s1,s2)
	fileOutput = open("output.txt", "w")
	outputString = str(r[2]) + '\n' + r[0] +'\n'+ r[1] +'\n'+str(format((time.time() - stime) * 1000,".6f")+'\n'+str(float(psutil.Process(os.getpid()).memory_info().rss /1024))+'\n')
	fileOutput.write(outputString)
	fileOutput.close()





