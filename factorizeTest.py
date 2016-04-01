"""
Author: Eric Knapik
Advanced Linear Algebra
QR Factorization in Python
"""


from numpy.linalg import inv
import numpy as np
import math
import os


"""
QR factorization implemented from:
https://en.wikipedia.org/wiki/QR_decomposition
http://www.math.ucla.edu/~yanovsky/Teaching/Math151B/handouts/GramSchmidt.pdf
"""
def qrFactorize(inputMat):

	matSize = inputMat.shape[0]
	# Initialize matricies
	matA = inputMat # set this so that I do not destroy the input matrix
	gn = np.identity(matSize)
	q = np.identity(matSize)

	step = 1.0
	val1 = 0.0
	val2 = 0.0
	cosX = 0.0
	sinX = 0.0

	i = 0
	while i < matSize:
		j = matSize - 1
		while j > i:
			value1 = matA[j-1, i]
			value2 = matA[j, i]
			cosX = value1 / math.sqrt((value1*value1)+(value2*value2))
			sinX = -value2 / math.sqrt((value1*value1)+(value2*value2))

			gn[j, j] = cosX
			gn[j, j-1] = sinX
			gn[j-1, j] = -sinX
			gn[j-1, j-1] = cosX

			#print("G" + str(step) + ":")
			#print(gn)
			matA = gn * matA

			#print("A" + str(step) + ":")
			#print(matA)
			q = q * gn

			gn = np.identity(matSize)

			step += 1
			j -= 1
		i += 1

	q = inputMat * inv(matA)
	return q, matA



def main():

	option = -1
	while (option < 0) or (option > 2):
		option = int(input("Options: \n(0) Quit\n(1) Use Hilbert Matrix\n(2) Enter your own nxn\n> "))
		try:
			os.system('cls')
			os.system('clear')
		except:
			pass

	if option == 0:
		return

	x = -1
	while x < 0:
		x = int(input("What is the size of the matrix: "))
	inputMat = np.mat(np.zeros((x,x)))

	if option == 1:
		# create a default matrix to use
		for i in range(0,x):
			for j in range(0,x):
				inputMat[i, j] = float(1.0 / ((i + 1.0) + (j + 1.0) - 1.0))
	else:
		for i in range(0,x):
			for j in range(0,x):
				inputMat[i, j] = float(input("Enter the value for (" + str(i) + ", " + str(j) +"): "))

	# All input entered now run the QR Factorization
	try:
		os.system('cls')
		os.system('clear')
	except:
		pass
	print("Input Matrix: ")
	print(inputMat)

	Q, R = qrFactorize(inputMat)
	print("Q: ")
	print(Q)
	print("R: ")
	print(R)
	print("Answer from QR multiplication: ")
	print(Q * R)

main()

