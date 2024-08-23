import numpy as np

class LocalAlignment:
	def __init__(self,gap,match,mismatch):
		self.gap = gap
		self.match = match
		self.mismatch = mismatch

	def misMatchChar(self,x,y):
		if x != y:
			return self.mismatch
		else:
			return self.match
	
	def getScore(self,xSeq, ySeq):
		score = 0
		for i in range (0, len(xSeq)):
			if xSeq[i] == ySeq[i]:
				score = score + self.match
			elif xSeq[i] == '-' or ySeq[i] == '-':
				score = score + self.gap
			elif xSeq[i] != ySeq[i]:
				score = score + self.mismatch
		return score

	def localAlign(self,x,y):
		# Create an initial matrix of zeros, its len(SeqA) x len(SeqB)
		matrix = np.zeros((len(x)+1, len(y)+1))
		traceback = np.zeros((len(x)+1, len(y)+1),dtype=str)
		traceback[0,0] = "f"
		best = 0
		optLoc = (0,0)

        # Fill in the matrix with alignment scores and obtain the best score and position
		for i in range(1,len(x)+1):
			for j in range(1,len(y)+1):
				left = matrix[i][j-1] + self.gap
				up = matrix[i-1][j] + self.gap
				diag = matrix[i-1][j-1] + (self.match if x[i-1] == y[j-1] else self.mismatch)
				matrix[i][j] = max(left,up,diag,0)

				if matrix[i][j] == left:
					traceback[i][j] = 'l'
				elif matrix[i][j] == up:
					traceback[i][j] = 'u'
				elif matrix[i][j] == diag:
					traceback[i][j] = 'd'
				else:
					traceback[i][j] = 0
		
				if matrix[i][j] >= best:
					best = matrix[i][j]
					optLoc = (i,j)
		
		return best, optLoc, matrix, traceback
	
	def getSequence(self,x,y):
		bestScore,optLoc,matrix,traceback= self.localAlign(x,y)
		# Obtaining the locally aligned sequence using matrix
		xSeq = []
		ySeq = []
		i = optLoc[0]
		j = optLoc[1]

		while(i > 0 or j > 0):

			diagonal = matrix[i-1][j-1]
			up = matrix[i-1][j]
			left = matrix[i][j-1]
	
			if traceback[i][j] == 'd':
				# Diag is scored when x[i-1] == y[j-1]
				xSeq.append(x[i-1])
				ySeq.append(y[j-1])
				i = i-1
				j = j-1
			elif traceback[i][j] == 'l':
				# Left holds true when '-' is added from x string and y[j-1] from y string
				xSeq.append('-')
				ySeq.append(y[j-1])
				j = j-1
			elif traceback[i][j] == 'u':
				# Up holds true when '-' is added from y string and x[j-1] from x string
				xSeq.append(x[i-1])
				ySeq.append('-')
				i = i-1
			elif traceback[i][j] == 'f':
				break
			else:
				# Break condition when diag score is the maximum
				break

		return int(bestScore), xSeq[::-1], ySeq[::-1]


# p1 = LocalAlignment(-3, 4, -1)
# x = "CTATTGACGTAACAT"
# y = "CTATTGAACAT"
# bestScore,xSeq, ySeq = p1.getSequence(x,y)

# print('The sequence obtained via traceback is: ')
# print(*xSeq)
# print(*ySeq)

# print('The best score obtained is: '+str(bestScore))