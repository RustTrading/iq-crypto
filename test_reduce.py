import iqlib, unittest, math
from random import randint
import cypari2 
import hashlib

class CryptoContext:

	def __init__(self, D):
		self.Discriminant_ = D
		self.BasePoint_ = iqlib.randomIdealFromPi(D)
		return

	public_keysB_ = []
	public_keysA_ = []
	private_keysB_ = []
	Discriminant_ = None
	BasePoint_ = None

def poly26(coeff):
	a,b,c = coeff
	return a * 26 ** 2 + b * 26 + c

def make_public_keysB(key_length, message_length, context):

	D = context.Discriminant_
	ideal = context.BasePoint_

	for _ in range (message_length):
		n = randint(1, 2 ** key_length)
		context.private_keysB_.append(n)
		ideal_public_keyB = iqlib.fast_reduce(iqlib.powern(ideal, n, D), D)
		context.public_keysB_.append(ideal_public_keyB)
	return 


def encrypt_messageA(message, key_length, context):

	ideal = context.BasePoint_
	D = context.Discriminant_

	message_num = [ ord(letter) - ord('A') for letter in message ]
	message_num_poly = []
		
	idx = 0
	r = len(message_num) % 3
		
	if r!=0:
		for _ in range(r, 3):
			message_num.append(25)

	public_keys = []		
	encrypted_message = []

	for k in range (len(message_num)//3):
		letter_tuple = (message_num[idx], message_num[idx + 1], message_num[idx + 2])
		poly_letter = poly26(letter_tuple)
		message_num_poly.append(poly_letter)
		n = randint(1, 2 ** key_length)
		idx = idx + 3
		ideal_public_keyAB = iqlib.fast_reduce(iqlib.powern(context.public_keysB_[k], n, D), D)

		context.public_keysA_.append(iqlib.fast_reduce(iqlib.powern(ideal, n, D), D))
		encrypted_message.append(poly_letter + ideal_public_keyAB[0])

	return encrypted_message


def decrypt_message(encrypted_message, context, key_length):

	idx = 0
	D = context.Discriminant_
	decrypted_message = []
	for pKeyB in context.private_keysB_:
		keyDecryptor = iqlib.fast_reduce(iqlib.powern(context.public_keysA_[idx], pKeyB, D), D)
		decrypted_message.append(encrypted_message[idx] - keyDecryptor[0])
		idx = idx + 1

	return decrypted_message


def recover_letter(letter_num):

	message = []
	letter_rest = letter_num

	res = letter_rest % 26
	message.append(chr(res + ord('A')))

	letter_rest = letter_num// 26

	for degree in range(1, 3):
		res = letter_rest % 26
		message.append('' + chr(res + ord('A')))
		letter_rest = letter_rest// 26

	message.reverse()
	return ''.join(message)


def recover_message(letters_encrypted):
	return [ recover_letter(l) for l in letters_encrypted ]


class TestReduceMethods(unittest.TestCase):
	#def test_egsd(self):
	#	print (egcd(96, 18))
	def test_identity(self):
		D = -55581980884703
		for a in range(1, int(math.sqrt(-D/3))):
			b = randint(1, a)
			if (b**2 - D) % (4* a) > 0 :
				continue
			i1 = (a , b)
			i2 = (a, -b)
			self.assertEqual(iqlib.fast_reduce(iqlib.multip(i1,i2,D), D) == (1,1), True)
		return	

	def test_multip0(self):
		D = -55581980884703
		i1 = (199, 2 * 190 + 1)
		i2 = (39601, 33415)
		i3 = (7880599, -5431523)


#		for i  in range(1, 5):
#			print(powern(i1, i, D))
		self.assertEqual(iqlib.fast_reduce(iqlib.multip(i1, i3, D), D) == iqlib.fast_reduce(iqlib.powern(i2, 2, D),D), True)
		return

	def test_mult_random(self):
		D = randint(10000, 1000000000000)
		pari = cypari2.Pari()
		
		dfactors = pari.factor(D)

		m = 5
		for d in dfactors[0]:
			if d > m and d % 4 == 3:
				m = int(d)

		print('Discriminant: ', -m)
		
		ideallist = []				
		for _ in range(1000):
			b = randint(1, 1000)
		
			if b % 2 == 0:
				b = b - 1
	
			A = ( b * b + m )//4		
			factorsa = pari.factor(A)
	
			for a in factorsa[0]:
				if a >= b:
					ideallist.append((int(a),b, A//int(a)))

		for _ in range(10):				
			idx1 = randint(0, len(ideallist))
			idx2 = randint(0, len(ideallist))

			form1 = pari.Qfb(*ideallist[idx1])
			form2 = pari.Qfb(*ideallist[idx2])

			form3 = form1 * form2

			a, b, c = form3

			id1 = (ideallist[idx1][0], ideallist[idx1][1])
			id2 = (ideallist[idx2][0], ideallist[idx2][1])

			self.assertEqual(iqlib.fast_reduce(iqlib.multip(id1, id2, -m), -m) == (a,b), True)
		return


	def test_multip1(self):
		a = 199
		b = 2 * 190 + 1
		D = -55581980884703

		ideal_reduced1 = iqlib.fast_reduce(iqlib.powern((a,b), 3, D), D)
		ideal_reduced2 = iqlib.fast_reduce(iqlib.powern((a,b), 7, D), D)

		ideal_reduced3 = iqlib.multip(ideal_reduced1, ideal_reduced2, D)
		ideal_reduced4 = iqlib.multip(ideal_reduced2, ideal_reduced1, D)
		
		self.assertEqual(ideal_reduced3 == ideal_reduced4, True)
		self.assertEqual(iqlib.check_if_reduced(ideal_reduced2, D), True)
		return

	def test_multip2(self):
		D = -143
		print(iqlib.multip((4,1), (8, -7), D))
		return

	def test_multip3(self):
		ideal = (2, 1)
		D = -143
		ideal_reduced = iqlib.reduce_bw(iqlib.powern(ideal, 3, D), D)
		print('1) Reduced:' , ideal_reduced)
		self.assertEqual(iqlib.check_if_reduced(ideal_reduced, D), True)
		return
	
	def test_ideal(self):
		ideal = (32, 43)
		D = -71
		self.assertEqual(iqlib.check_if_reduced(ideal, D), False)
		return
	
	def test_power(self):
		ideal = (4, 3)
		D = -71
		ideal_power = iqlib.reduce_bw(iqlib.powern(ideal, 50, D), D)
		print('2) Reduced:' ,ideal_power)
		self.assertEqual(iqlib.check_if_reduced(ideal_power, D), True)
		return

	def test_reduce1(self):
		D = -143
		power_ideal = iqlib.powern((2, 1), 10, D)
		ideal_reduced = iqlib.reduce_bw(power_ideal, D)
		a, b = power_ideal
		c = (b * b - D)//(4 * a)
		print ('3) Reduced:', ideal_reduced)
		self.assertEqual(iqlib.check_if_reduced(ideal_reduced, D), True)
		power_form_reduced = iqlib.form_reduce((a,b,c))
		self.assertEqual(ideal_reduced == (power_form_reduced[0], power_form_reduced[1]), True)
		#print(reduce(ideal, delta))
		return

	def test_reduce2(self):
		ideal = (32, 43)
		D = -71
		ideal_reduced= iqlib.reduce_bw(ideal, D)
		print('4) Reduced: ', ideal_reduced)
		self.assertEqual(iqlib.check_if_reduced(ideal_reduced, D), True)
		return

	def test_form_reduce(self):
		form = (195751, 37615, 1807)
		a, b, c = form
		D = b * b - 4 * a * c
		form_reduced = iqlib.form_reduce(form)
		print ('5) Reduced: ', form_reduced)
		self.assertEqual(iqlib.check_form_reduced(form_reduced), True)
		ideal_reduced = iqlib.fast_reduce((a,b), D)
		self.assertEqual(ideal_reduced == (form_reduced[0], form_reduced[1]), True)
		return

	def test_rho(self):
		form = (195751, 37615, 1807)
		print ('Rho test', iqlib.rho_form(form))
		return

	def test_classgroup(self):
		c1 = (2, 1)	
		D = -143
		print(iqlib.fast_reduce(iqlib.powern(c1, 6, D),D))
		return

	def test_random_ideal(self):
		D = -32127240703
		order1, order2, order3 = 1307 , 11 , 7
		orderClass = 100639
	
		random1 = iqlib.randomIdealFromPi(D)
		random2 = iqlib.randomIdealFromPi(D)
		random3 = iqlib.randomIdealFromPi(D)

		random4 = iqlib.multip(random1, random2, D)
		random4 = iqlib.multip(random4, random3, D)
    
		pari = cypari2.Pari()

		print ('Random Ideals: ', random1, random2, random3)
		print( 'ramdom1 a factors: ', pari.factor(random1[0]))
		print( 'ramdom2: a facors ', pari.factor(random2[0]))
		print( 'ramdom3: a factors', pari.factor(random3[0]))

		print ('Order1: ', random1, iqlib.terrAlgo(random1, D))
		print ('Order2: ', random2, iqlib.terrAlgo(random2, D))
		print ('Order3: ', random3, iqlib.terrAlgo(random3, D))
		print ('Order random1*random2*random3: ', random3, iqlib.terrAlgo(random4, D))

		print(iqlib.fast_reduce(iqlib.powern(random1, orderClass, D),D))
		print(iqlib.fast_reduce(iqlib.powern(random2, orderClass, D),D))
		print(iqlib.fast_reduce(iqlib.powern(random3, orderClass, D),D))

		print(iqlib.fast_reduce(iqlib.powern(random1, order1, D),D))
		print(iqlib.fast_reduce(iqlib.powern(random2, order2, D),D))
		print(iqlib.fast_reduce(iqlib.powern(random3, order3, D),D))
		return
	
	def  test_encrypt_decrypt_prime(self):

		D = -iqlib.getRandomPrime(2**1024, 2**2048)
		print ('Discriminant: ', D, D % 4 )
		context = CryptoContext(D)
		ideal = context.BasePoint_
		print ('Base Point Ideal: ', ideal)

		message = 'Happy Thanksgiving and Merry Christmas and The New Year and Other Celebrations'
		message = message.replace(' ', '')
		message =  message.upper()

		make_public_keysB(100, (len(message) + 2)//3, context)
		#print (context)
		print ('encrypting....')
		encrypted_message = encrypt_messageA(message, 100, context)
	
		print ('Encryption: ', encrypted_message)
		print ('Public KeysA: ', context.public_keysA_)
		print ('Public KeysB: ', context.public_keysB_)
		print ('Decryption: ', ''.join(recover_message(decrypt_message(encrypted_message, context, 100))))
		return
	
	'''	
	def test_sign_message(self):
		
		D = iqlib.getRandomPrime(2**128, 2**256)

		print ('Discriminant: ', -D)		

		key = str(randint(1, 2**128))

		keyhash = hashlib.sha256(key.encode('utf-8'))
		ideal_power = int.from_bytes(keyhash.digest(), 'big')
		
		print('PrivateKey: ', key)
		print ('KeyHash', keyhash.digest())
		print ( 'ideal: power: ', ideal_power)
	'''

	def test_classnumber(self):

		deltas = [ -3, -4, -7, -8, -11, -15, -19, -20, -23, -24, -28, -31, -163, -79 ]
		classDim = [ 1, 1, 1, 1, 1, 2, 1, 2, 3, 2, 1, 3, 1, 5 ]

		for d in enumerate(deltas):
			self.assertEqual(iqlib.classNumber(d[1]) == classDim[d[0]], True)
		return
			
	'''
	def test_encrypt_decrypt(self):
		D = -55581980884703
		ideal = (199, 2 * 190 + 1)

		message = 'Happy Thanksgiving and Merry Christmas and The New Year and Other Celebrations'
		message = message.replace(' ', '')
		message =  message.upper()

		context = CryptoContext()
		make_public_keysB(ideal, D, 100, (len(message) + 2)//3 , context)
		encrypted_message = encrypt_messageA(message, ideal, D, 100, context)
	
		print ('Encryption: ', encrypted_message)
		print ('Public KeysA: ', context.public_keysA_)
		print ('Public KeysB: ', context.public_keysB_)
		print ('Decryption: ', ''.join(recover_message(decrypt_message(encrypted_message, context, 100))))
	'''	
	'''	m, counter = 0, 0 

		while counter < 10:
			random = iqlib.randomIdealFromPi(D)
			n = iqlib.terrAlgo(random, D)
			if n > m :
				m = n
			elif m % n == 0:
				counter = counter + 1
				if n < m:
					print ('Smaller Order: {0}, Class Number update: {1}'.format(n, m))	
			else:
				d = iqlib.gcd(m,n)
				m = (n/d) * m

		for counterCheck in range(10):
			random = iqlib.randomIdeal(D)
			ideal = iqlib.fast_reduce(iqlib.powern(random, m, D),D)
			self.assertEqual(ideal == (1,1), True)
					
		print ('Class Number : ', m)	
		'''	
if __name__ == '__main__':
    unittest.main()
