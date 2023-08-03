'''
	pi number 10^9 digits: https://stuff.mit.edu/afs/sipb/contrib/pi/pi-billion.txt

strong prime:  9223026960492503926560447883691216823891159206945879338569055
6112885081234276554752432662947430471254857392209173424738670836432261695415
2696342713007777237076609751362536974170509151486476053780070525961519300239
2925974322799881779279826923181567660057189570005135889154307919170273097771
2806257105999935839 offset:  670644371 bitlength:  1024
residue: p % 4 3
Pari Prime:  1
strong prime:  90564999926208252939531685806691926244923520267874489353411975
97549226344029015963116572165125063925531210008256518660363155248570415696296
70501565179707678173469574095566676758254848516876040457986939777368793971602
09246823165439340225219102360915477043504922953761823150898696611455181533811
737441737521737 offset:  849714871 bitlength:  1024
residue: p % 4 1

BasePoint Ideal:  (90564999926208252939531685806691926244923520267874489353411
	9759754922634402901596311657216512506392553121000825651866036315524857041
	5696296705015651797076781734695740955666767582548485168760404579869397773
	68793971602092468231654393402252191023609154770435049229537618231508986966
	11455181533811737441737521737, 
	45428565791202071817923471243695412481527232151647722274787277130389257972
	40158416955068658241931299523354338865723744325780819954618564816558689095
	60850550754166641775736230978678617193488796549511824695179351955912276844
	80795456742795909542847395767928742798330764194373793483203350254654337901
	021102980529)
'''

import cypari2
from random import randint, seed
from math import sqrt, log10, log2
import iqlib

def egcd(a, b):
	x,y, u,v = 0,1, 1,0
	while a != 0:
		q, r = b//a, b % a
		m, n = x-u*q, y-v*q
		b,a, x,y, u,v = a,r, u,v, m,n
	gcd_val = b
	return gcd_val, x, y


def modInv(x, p):
	_ , xInv, _ = egcd(x, p)
	if xInv < 0:
		xInv = p + xInv
	return xInv	


def mult_mont(aMont, bMont, n, r, nprime):

	t = aMont * bMont
	m = t * nprime % r
	u = (t + m * n)//r
	#print (aMont, ',', bMont, ',', t,',', m, ',', u)
	return u % n


def exp_mont(M, e, n):
	k = M.bit_length() + 1
	r = 2 ** k
#	print('r: ', r)
	rInv = modInv(r, n)
#	print('r^-1: ', rInv)

	nprime =(r * rInv -1)//n

#	print ('nprime: ', nprime)

	MMont = M * r % n
	xMont = r % n

	i = k - 1
	while i >=0 :
#		print('Step 5 before: ', xMont)
		xMont = mult_mont(xMont, xMont, n, r, nprime)
#		print('Step 5: ', xMont)
		if e & 2**i > 0:
#			print('Step 6 before: ', MMont, '*', xMont)
			xMont = mult_mont(MMont, xMont, n, r, nprime)
#			print('Step 6: ', xMont)
		i = i - 1
	x = mult_mont(xMont , 1, n, r, nprime)	
	return x


def miller_rabin(n, r = 25):
	k = n - 1
	s = 0

	while True:
		if k > 0 and k % 2 == 0 :
			s = s + 1
			k = k//2
		else:
			break
	#print (n, ' = 1 + 2 ^ {0} * {1}'.format(s, k))		

	witnessLoop = False

	for _ in range(r):

		a = randint(2, n-2)
		x = exp_mont(a, k, n)

		if x == 1 or x == n - 1:
			continue

		for _ in range(s - 1):

			x = (x % n)
			x = x * x % n

			if x == 1 :
				return False
			
			if x == n - 1:
				witnessLoop = True
				break
		
		if not witnessLoop:
			return False

	return True			


def generate_a(fixLength):
	#g = open('pi-billion.txt', 'r')
	g = open('hash_number.txt', 'r')
	v = g.read()
	s = len(v)
	
	pari = cypari2.Pari()
	strongPrime =False
	
	primeLength = int( (fixLength - 1) * log10(2))
	variationLength = 2
	maxLength = primeLength + variationLength
	
	seed()

	while True:

		i = randint(2, s - maxLength)
		j = i + primeLength

		for k in range(variationLength) :

			p = int(v[i : j + k])

			if  p.bit_length() == fixLength :
				break
				
		if  p.bit_length() != fixLength :
				continue		

		if miller_rabin(p):
			
			print('strong prime: ', p, 'offset: ', i, 'bitlength: ', p.bit_length())
			print('residue: p % 4', p % 4)
			print ('Pari Prime: ', pari.isprime(p))

			return p

	return None


def generate_d(fixLength):

	g = open('pi-billion.txt', 'r')
	v = g.read()
	s = len(v)
	p = 2

	pari = cypari2.Pari()
	strongPrime =False
	
	primeLength = int((fixLength - 1) * log10(2))
	maxLength = primeLength
	variationLength = 2

	seed()

	while True:

		i = randint(2, s - maxLength)
		j = i + primeLength
		
		tail = int(v[j - 2 : j])

		if tail % 4 != 3 :
			continue

		for k in range(variationLength) :

			p = int(v[i - k : j])

			if  p.bit_length() == fixLength :
				break

#		while j - i <  3 * primeLength // 2 :

#		p = int(v[i : j])
		
		if  p.bit_length() != fixLength :
				continue		
	
		if miller_rabin(p):
			
			print('strong prime: ', p, 'offset: ', i, 'bitlength: ', p.bit_length())
			print('residue: p % 4', p % 4)
			print ('Pari Prime: ', pari.isprime(p))
				
			return p	

		#	if i  > 2 :
		#		i = i - 1
		#	else :
		#		break	

	return None			
		

def generate_commitment(n):
	
	d = generate_d(1024)
	primelist = []

	while True:
		
		prime = generate_a(1024)
		
		if iqlib.LegendreSymbol(-d, prime) + 1 <= 0:
			continue

		if len(primelist) < n :
			primelist.append(prime)
		else:
			return (-d, primelist)

	return None	


			