import numpy as np

class GlobalAlignment:
	def __init__(self,gap,match,mismatch):
		self.gap = gap
		self.match = match
		self.mismatch = mismatch
	
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

	def getMatrix(self,x,y):
		# Create an initial matrix of zeros, its len(SeqA) x len(SeqB)
		sizeX = len(x)
		sizeY = len(y)
		matrix = np.zeros((sizeX+1, sizeY+1))
		traceback = np.zeros((sizeX+1, sizeY+1),dtype=str)

		# Initializing the first cell in row and column with done--> first cell (o)
		traceback[0][0] = 'f'
		# Initializing the first row and first column with the gap values in matrix and up(u) or left(l) values in traceback matrix
		for j in range(1,sizeY+1):
			matrix[0][j] = j*self.gap
			traceback[0][j] = 'l'
		for i in range(1,sizeX+1):
			matrix[i][0] = i*self.gap
			traceback[i][0] = 'u'
		
		return matrix,traceback
	
	def globalAlign(self,x,y):
		# Fill in the matrix with alignment scores
		matrix,traceBack = self.getMatrix(x,y)
		
		for i in range(1,len(x)+1):
			for j in range(1,len(y)+1):
				left = matrix[i][j-1] + self.gap
				up = matrix[i-1][j] + self.gap
				diag = matrix[i-1][j-1] + (self.match if x[i-1] == y[j-1] else self.mismatch)
				matrix[i][j] = max(left,up,diag)
				if matrix[i][j] == left:
					traceBack[i][j] = 'l'
				elif matrix[i][j] == up:
					traceBack[i][j] = 'u'
				else:
					traceBack[i][j] = 'd'
		return matrix,traceBack
	
	def getAlignedSequences(self,x,y):
		matrix,traceBack = self.globalAlign(x,y)
		# Obtain x and y globally aligned sequence arrays using the bottom-up approach
		xSeq, ySeq = [], []
		i, j = len(x), len(y)
		while(i > 0 or j > 0):
			if traceBack[i][j] == 'd':
				# Diag is scored when x[i-1] == y[j-1]
				xSeq.append(x[i-1])
				ySeq.append(y[j-1])
				i = i-1
				j = j-1
			elif traceBack[i][j] == 'l':
				# Left holds true when '-' is added from x string and y[j-1] from y string
				xSeq.append('-')
				ySeq.append(y[j-1])
				j = j-1
			elif traceBack[i][j] == 'u':
				# Up holds true when '-' is added from y string and x[j-1] from x string
				xSeq.append(x[i-1])
				ySeq.append('-')
				i = i-1
			elif traceBack[i][j] == 'f':
				# Break condition when we reach the [0,0] cell of traceback matrix
				break
		score = self.getScore(xSeq, ySeq)
		return score, xSeq[::-1], ySeq[::-1]
	

# p1 = GlobalAlignment(-3, 4, -1)
# # x = "aaac"
# # y = "agc"
# x = "CTATTGACGTAACAT"
# y = "CTATTGAACAT"
# print(x)
# score, xSeq, ySeq= p1.getAlignedSequences(x,y)
# print('The globally aligned sequences are:')
# print(*xSeq)
# print(*ySeq)
# print("score = ", score)