#   Assume everywhere that D < 0 and D = 1 (mod 4)
#
#   Definition:
#	A positive binary quadratic form ax^2 + bxy + c^2  is normal if -a < b <= a
#	
#	here c = (b^2 - D)/4a
#   Definition:
#	The positive definite form (a, b, c) is called reduced if it is normal, a <= c, and if b >= 0 for a = c
#	We assume that w = (1 + sqrt(-D))/2, that is r = 2 in w = ( r - 1 + sqrt(-D)) / r

from math import frexp, sqrt, log10, log2
from random import randint
from mpmath import mp
import decimal
import cypari2

def check_if_reduced(ideal, D):
	a, b = ideal
	if ( b * b - D ) % (4 * a) != 0:
		return False
	c = ( b * b - D )//(4 * a)
	return check_form_reduced((a, b, c))


def check_form_reduced(form):
	a, b, c = form

	if abs(b) > a or b == -a:
		return False

	if a > c or (a == c and b < 0):
		return False

	return True	


def egcd(a, b):
	x,y, u,v = 0,1, 1,0
	while a != 0:
		q, r = b//a, b % a
		m, n = x-u*q, y-v*q
		b,a, x,y, u,v = a,r, u,v, m,n
	gcd_val = b
	return gcd_val, x, y


def gcd(u, v):
	while v:
		u, v = v, u % v
	return abs(u)


def multip(ideal1, ideal2, D):
	a1, b1 = ideal1
	a2, b2 = ideal2

	d1, x1, y1 = egcd(a1, a2)
	a3 = a1 * a2
	d2, x2, y2 = egcd((b1 + b2)//2, d1)

	if d1 == d2 :
		x2, y2 = 0 ,1

	a3 = a3//(d2 * d2)
	b3 = (y2 * x1 *( b2 - b1 ) + x2 * (D - b1 * b1)//(2 * a1)) % (2 * a2 // d2)
	b3 = b1 + b3 * (a1//d2) % (2 * a3)
	if abs(b3) > a3:
		b3 = b3 - 2 * a3 * nearest_int(b3, 2 * a3)
	return (a3, b3)


def nearest_int(P, Q):
	q, R = P//Q, P % Q
	if 2 * R >= Q:
		return q + 1
	return q	


def reduce_bw(ideal, D):
	Q0, P0 = ideal
	Q0 = 2 * Q0
	while True:
		q0 = nearest_int(P0, Q0)
		P1 = q0 * Q0 - P0
		Q1 = (P1 * P1 - D)//Q0
		if Q1 >= Q0:
			P0 = P0 % Q0
			if P0 >= Q0//2:
				P0 = P0 - Q0

			if Q0 //2 == -P0 or (Q0* Q0 == (P0 * P0 - D) and P0 < 0 ):
				P0 = -P0	
			return (Q0//2, P0)

		Q0 = Q1
		P0 = P1

	return (Q0//2, P0)


def rho_form(form):
	a, b, c = form
	s = nearest_int(b, 2 * c)
	return (c, -b + 2 * s * c, c * s * s - b * s + a)


def fast_reduce(ideal, D):
	return reduce_bw(ideal,D)
	a, b = ideal
	c = ( b * b - D )//(4 * a)
	a, b, _ = form_reduce((a, b, c))
	return (a, b)


def form_reduce(form):
	a, b, c = form
	form_normalized = form

	while not check_form_reduced(form_normalized):
		form_normalized = rho_form(form_normalized)

	return form_normalized


def power2(ideal, D):
	a1, b1 = ideal
	a3 = a1 * a1
	d2, x2, y2 = egcd(b1, a1)

	if a1 == d2 :
		x2, y2 = 0 ,1

	a3 = a3//(d2 * d2)
	b3 = (x2 * (D - b1 * b1)//(2 * a1)) % (2 * a1 // d2)
	b3 = b1 + b3 * (a1//d2) % (2 * a3)
	if abs(b3) > a3:
		b3 = b3 - 2 * a3 * nearest_int(b3, 2 * a3)
	return (a3, b3)


def LegendreSymbol(D,p):

	if p == 2:
		if D % 2 == 0:
			return 0
		else:
			return (-1) ** ((D ** 2 - 1)//8)

	ls = exp(D, (p-1)//2,  p)
	if ls == p - 1:
		ls = -1
	return ls	


def modInv(x, p):
	_ , xInv, _ = egcd(x, p)
	if xInv < 0:
		xInv = p + xInv
	return xInv


def numberOfPrimeForms(D, p):
	return LegendreSymbol(D, p) + 1


def sqrtModP(r, p):
	m = p-1
	t = 0
	while m % 2 == 0:
		t = t + 1
		m = m // 2

	for _ in range(p):
		c = randint(1, p)
		if LegendreSymbol(c, p) == 1 :
			continue
		rm = exp(r, m, p)
		cm = modInv(exp(c, m, p), p)
		e, i = 0, 0	
		while rm !=1 :
			i = i +1
			cm = cm**2
			if exp(rm, 2 ** (t-i-1), p)!=1:
				e = e + 2**i
				if rm != rm * cm % p :
					rm = rm * cm % p
				else:
					break	

		if rm !=1:
			continue

		a = modInv(exp(c, e, p), p)
		a = a * r			
		return exp(c, e//2,  p) * exp(a, (m+1)//2, p) % p

	return None	


def sqrtMod4P(D, p):
	if p == 2:
		if D % 2 == 0:
			return 2 * ( D//4 % 2)
		else:
			return 1
	else:
		r =	sqrtModP(D, p)
		if r == None:
			return None

		if (r - D) % 2 == 0:
			return r
		else:
			return p-r


def primeFormNumber(D, p):
	return LegendreSymbol(D,p) + 1


def primeForm(D, p):

	if primeFormNumber(D, p) == 0:
		return None

	b = sqrtMod4P(D,p)

	if b == None:
		return None

	return (p, b, (b**2 - D)//(4 * p))


def classNumber(D):
	h = 0
	b = D % 2
	while b <= sqrt(-D/3):
		A = (b * b -D)//4
		a = max(1,b)

		while a <= sqrt(A):
			if A % a == 0:
				c = A//a
				if gcd(a, gcd(b ,c)) == 1:
					if b == 0 or a == b or a == c :
						h = h + 1
					else:
						h = h + 2		
			a = a + 1
		b = b + 2	
	return h		


def getRandomPrime(start, end):
	key = 0
	pari = cypari2.Pari()
	while True:
		key = randint(start, end)
		prime = int(pari.nextprime(key))
		if  prime % 4 == 3 :
			return prime


def randomIdeal(D):

	dec = int(decimal.Decimal(-D//3).sqrt())
	p = getRandomPrime(1, dec)

	while primeFormNumber(D, p) == 0 :
		p = getRandomPrime(1, dec)

	pform = None
	counter = 0

	while pform == None:
			pform = primeForm(D, p)
			counter = counter + 1
			if counter > 1000:
				counter = 0
				p = getRandomPrime(1, dec)
	
	if (pform[1] * pform[1] - D) % (4 * p) > 0:
		print ('Wrong Ideal Generated: ', pform)
		raise ValueError()

	return (pform[0], pform[1])		

def getEssentialIdeals(D):
	B = 6 * log2(-D)**2
	pari = cypari2.Pari()
	p = 2
	while p < B:
		print(primeForm(D, int(p)))
		print ('prime: ', p)
		p =  pari.nextprime(p+1)

	return	
	

def randomIdealFromPi(D):

	mp.dps= 1000
	a = b = None
	key_degree = int(log10(-D//3)//2)//10000
	
	dec = int(decimal.Decimal(-D//3).sqrt())
	pari = cypari2.Pari()

	while True:
		st = randint(1, 990)
		b = int(10 ** key_degree * (mp.pi * 10 ** st  - int(mp.pi * 10 ** st)))

		while b > dec :
			st = randint(1, 990)
			b = int(10 ** key_degree * (mp.pi * 10 ** st  - int(mp.pi * 10 ** st)))
		
		if b % 2 == 0:
			b = b -1


		A = (b * b - D)//4
		a == 1
		
		print('test', A, b)

		while a == 1 :
			a = randint(b, int(decimal.Decimal(A).sqrt()))
			a = gcd(a, A)
			print('test', a)

		c = A//a
		
		if gcd(a, gcd(b ,c)) == 1:
			return (a, b)
	
	return None


def exp(n, m, modp):

	if m == 0:
		return 1
	
	if m == 1:
		return n % modp

	rest = m
	num = rest.bit_length() - 1
	powerm = n
	powerm_prev = 1
	while num > 0:	
		for _ in range(num):
			powerm = powerm * powerm % modp

		powerm = powerm_prev * powerm % modp
		rest = rest - 2 ** num
		if rest == 0 or rest == 1:
			break	
		num = rest.bit_length()-1
		powerm_prev = powerm
		powerm = n

	if m % 2 == 1:
		powerm = n * powerm % modp

	return powerm


def powern(ideal, n, D):

	if n == 0:
		return (1 , 1)

	if n == 1:
		return ideal

	rest = n
	num = frexp(rest)[1]-1
	ideal_power = fast_reduce(ideal, D)
	ideal_prev = (1, 1)
	while num > 0:	
		for _ in range(num):
			ideal_power = fast_reduce(power2(ideal_power, D),D)

		ideal_power = fast_reduce(multip(ideal_prev, ideal_power, D),D)
		rest = rest - 2 ** num
		if rest == 0 or rest == 1:
			break	
		num = frexp(rest)[1]-1
		ideal_prev = ideal_power
		ideal_power = ideal

	if n % 2 == 1:
		ideal_power = fast_reduce(multip(ideal, ideal_power, D),D)

	return ideal_power


def terrAlgo(ideal, D):
	babySet = {}
	babySet[(1,1)] = 0
	e = 1
	babyElement = ideal
	giantElement = ideal
	while True:
		if giantElement in babySet :
			return e * (e + 1)//2 - babySet[giantElement]
		babySet[babyElement] = e
		babyElement = fast_reduce(multip(ideal, babyElement, D),D)
		e = e + 1
		giantElement = fast_reduce(multip(giantElement, babyElement, D),D)
		