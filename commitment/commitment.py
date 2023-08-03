import pideal, iqlib
import cypari2
import decimal
from math import log

def main():

	pari = cypari2.Pari()
	idealData = pideal.generate_commitment(8)

	print('Ideal data: ', idealData)

	D, plist = idealData

	ilist = []

	for p in ilist:
		a, b, _ = iqlib.primeForm(D, p)
		ilist.append((a,b))

	delta = decimal.Decimal(-D)
	
	n = int((delta.ln()/2 + delta.ln().ln() + 1) *  delta.sqrt())

	print('p order: ', pari.nextprime(n))
	print('g, h - Ideal: ', ilist)

	return

def make_pprod(glst, hlst, u, alst, blst, c):

	gprod = None
	hprod = None

	for i in range(len(glst)):
		gprod = iqlib.fast_reduce(iqlib.multip(iqlib.multip(iqlib.powern(glst[i], alst[i], D), gprod, D), D), D)

	for i in range(len(hlst)):
		hprod = iqlib.fast_reduce(iqlib.multip(iqlib.multip(iqlib.powern(hlst[i], blst[i], D), hprod, D), D), D)
	
	return iqlib.fast_reduce(iqlib.multip(iqlib.multip(gprod, hprod, D), iqlib.powern(u, c, D) ,D), D)


def commit(g, h, u, P, a, b):
	return


def open(a, b, g, h, u):
	
	D = 123
    p = 2

	n = len(a)

	if n == 1 :
		c = a * b
		return P == iqideal.multip(iqideal.multip(iqideal.powern(g[0], a, D), iqideal.powern(h[0], b, D), D), iqideal.powern(u, c, D))
	else:
		n = n//2
		cl = numpy.dot(a[:n], b[n:]) % p
		cr = numpy.dot(a[n:], b[:n]) % p
		L = 1
		R = 2

if __name__ == '__main__':
	main()
