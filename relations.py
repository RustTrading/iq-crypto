import iqlib
from math import log2
from random import randint
import cypari2
import itertools

'''
Calculate relations
D = -3299

'''

def main():
		#D = -18285
		#D = -3299
		D = -55581980884703
		#D = -32127240703
		#D = -143
		#D = -103387
		B = 6 * int(log2(-D)**2)

		pari = cypari2.Pari()
		p = 2
		
		print('ClassBound: ', B)
		idealList = []
		primeIdeals = {}

		while p < B :
			id = iqlib.primeForm(D,int(p))
			p = pari.nextprime(p+1)

			if id != None:
				a,b, _ = id
				idReduced = iqlib.fast_reduce((a,b), D)
				idealList.append(idReduced)


		reduced = set(idealList)
	
		primeReduced = []

		for f in reduced:
			if pari.isprime(f[0]):
				primeReduced.append(f)

		ideals_len = randint(1, len(primeReduced))

		print('Ideals: ', ideals_len, len(primeReduced))
		
		ideals = set()

		for i in range(ideals_len):
			ideals.add(primeReduced[randint(0, ideals_len-1)])

		relations = set()

		ideals = list(ideals)

		count_relations = 0

		threshHold = B//10

		for _ in range(100000):
			
			x = randint(0, B)
			isTrivial = (x == 0)
			ideal_current = iqlib.fast_reduce(iqlib.powern(ideals[0], x, D),D)

			relation = []
			relation.append(x)
			prod = None

			for ideal in ideals[1:]:
				x = randint(0, B)
				
				if isTrivial:
					isTrivial = (x == 0)

				relation.append(x)	
				prod = iqlib.fast_reduce(iqlib.multip(ideal_current, iqlib.powern(ideal, x, D),D),D)
				
				#print ('Prod: ', prod[0], threshHold)

			if 	prod == (1,1) and not isTrivial:
				count_relations = count_relations + 1
				print('relation: ', relation)

				if count_relations > 25:
					break

				relations.add(tuple(relation))
				#print ('Prod: ', prod[0], threshHold)

			elif prod[0] < threshHold and not isTrivial:
				
				factors = pari.factor(prod[0])

				#print('factors: ', factors)
				#print ('prime ideals: ', primeIdeals)

				count = -1
				ideal2 = None

				for f in factors[0]:
					count = count + 1
					if int(f) in primeIdeals:
						ideal2 = primeIdeals[int(f)]
					else:
						continue	
					
					power = factors[1][count]
					a,b = iqlib.powern(ideal2, count, D)
					prod = iqlib.fast_reduce(iqlib.multip(prod, (a, -b),D), D)

					#print('factorization: ', prod)

				if prod == (1,1):
					print('discovered relation: ', relation)
					relations.add(tuple(relation))
					count_relations = count_relations + 1


		print ('Relations', relations)
		v = list(itertools.chain(*relations))

		print (v)

		m = pari.matrix(len(relations), len(ideals), v)
		smith = pari.matsnf(m)


		for i in range(len(smith)):
			if smith[i] > 1:
				print ('generator vs order: ', ideals[i % len(ideals)], smith[i])
				print ('Check: ', iqlib.terrAlgo(ideals[i % len(ideals)], D))

				
if __name__ == '__main__':
	main()					