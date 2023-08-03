import numpy as np
from math import sqrt; from itertools import count, islice
#from iqcrypto.test_reduce import egcd
#from test_reduce import egcd
from primesieve import primes
from decimal import Decimal
import random
import math
import string
from datetime import datetime
import cypari2
from numpy import dot, zeros
from numpy.linalg import matrix_rank, norm
import time
import copy
import sympy as sym
class ideal:
    def __init__(self, a, b, D):
        self.a = int(a)
        self.b = int(b)
        self.D = D
        self.c = (self.b*self.b - self.D)/(4*self.a)
        if self.c == int(self.c):
            self.c = int(self.c)
        else:
            print("c = ", self.c, " is not ideal")

    def deltaF(self):
        if self.D % 4 == 1:
            return self.D
        else:
            return 4 * self.D

    def is_reduced(self):
        if (abs(self.b) <= self.a) and (self.b*self.b - self.deltaF() >= 4*self.a*self.a):
            return True
        else:
            return False
    def iprint(self, str):
        print(str, " a = ", self.a, " b = ", self.b, " c = ", self.c, "D = ", self.D, " reduced ", self.is_reduced(), " prime ", self.is_prime(), " inversible ", self.is_inversible(), " ideal ", self.is_ideal())

    def is_prime(self):
        if is_prime(self.a) or self.a == 1:
            if (legendre_symbol(self.D, self.a) != -1):
                #print("D != a", 0 <= self.b, self.b <= self.a, self.a, self.b)
                if (0 <= self.b and self.b <= self.a):
                    #print("|b|<= a", -self.D %(4*self.a))
                    if self.b* self.b %(4*self.a) == self.D %(4*self.a):
                        return True
        return False
    def norm(self):
        #return int(fast_reduce(self).a)
        return self.a

    def element_norm(self, x, y):
        a = self.a
        b = self.b
        D = self.D
        return (a*a*x*x + a*b*x*y + (b*b - D)*y*y/4)

    def F(self, x):
        #print(" x ", x, " a ", self.a, " b ", self.b, " c ", self.c)
        return self.a*x*x + self.b*x + self.c

    def f(self, x, y):
        return self.a*x*x + self.b*x*y + self.c*y*y

    def is_inversible(self):
        gcd_bc = egcd(self.b, self.c)
        gcd_abc = egcd(self.a, gcd_bc[0])
        if gcd_abc[0] == 1:
            return True
        else:
            return False

    def is_ideal(self):
        if self.a//1 == self.a and self.b//1 == self.b and self.c//1 == self.c:
            return True
        else:
            return False

    def conjugation(self):
        return ideal(self.a, -1 * self.b, self.D)

    def integer_mul(self, int):
        self.a = int * self.a
        self.b = int * self.b
        self.c = int * self.c


def r(D):
    if D.real%4 == 1:
        r = 2
    else:
        r = 1
    return r


def starting_ideal(D):
    while True:
        q = int(np.random.uniform(1, 1000))
        # q = 5
        a = 4 * q + 3
        if not is_prime(a):
            continue
        if (legendre_symbol(D, a) == 1):
            break
    if D % 4 != 1:
        prime_ideal = primeForm(D, a)
        bp = prime_ideal.b
        b = bp + 2*a*np.random.randint(-10000, 10000)
    else:
        power = int((a + 1) / 4)
        b = D % (r(D) * a)
        b_temp = 1
        for i in range(int(power)):
            b_temp *= b
            b_temp = b_temp % (r(D) * a)
        b = b_temp
        b = int((b - r(D) + 1) / r(D))
        b = 2 * b + 1  # convert from form (a, b + omega) to (a, (b+sqrt(D))/2) form
    return ideal(a, b, D)



def is_prime(n):
    return n > 1 and all(n%i for i in islice(count(2), int(sqrt(n)-1)))

def legendre_symbol(a, p):
    if p == 2:
        if a % 2 == 0:
            return 0
        else:
            return (-1)**((a*a - 1)/8)
    #print("D ", a, " a ", p)
    modp = np.mod(a, p)
    if (modp == 0):
        ls = 0
        return ls
    power = int((p - 1)/2)
    for i in range(1, power):
        modp = (modp*np.mod(a, p)) % p
    if modp == 1:
        ls = 1
    else:
        ls =-1
    return ls

def convert_massage(str):
    str = str.replace(' ', '')
    di = dict(zip(string.ascii_letters, [(ord(c) % 32)-1 for c in string.ascii_letters]))
    step = 3
    base = 26
    ibeg = 0
    iend = step
    if(len(str)%step == 0):
        size = len(str)//step
        dummysize = 0
    else:
        size = (len(str)//step)+1
        dummysize = size * step - len(str)

    for i in range(dummysize):
        str += "z"
    m = []
    for i in range(size):
        substr = str[ibeg:iend]
        M = 0
        for j, char in enumerate(substr):
            power = j + 1
            M += di[char]*pow(base, step - power)
        ibeg += step
        iend += step
        m.append(M)
    return m



def reduce(idealA):
    if idealA.is_reduced():
        return idealA
    Q0 = idealA.a
    P0 = idealA.b
    D = idealA.D
    deltaF = idealA.deltaF()
    Q0 = 2 * Q0
    while True:
        q0 = (2 * P0 + Q0) // (2 * Q0)
        P1 = q0 * Q0 - P0
        Q1 = (P1 * P1 - deltaF) // Q0
        if Q1 >= Q0:
            P0 = P0 % Q0
            if P0 >= Q0 // 2:
                P0 = P0 - Q0
            return ideal(Q0 // 2, P0, D)
        Q0 = Q1
        P0 = P1
    return ideal(Q0 // 2, P0, D)


def iexp(idealA, n):
    ideal0 = idealA
    #print("kkkkscdfjdf")
    D = idealA.D
    if n == 0:
        return I(D)
    if n < 0:
        ideal0.b = - ideal0.b
    N = n
    imax = 0
    while int(N / 2) != 0:
        imax += 1
        N = N / 2
    N = abs(n)
    p = [True] * (imax+1)
    ideals = [ideal0] * (imax+1)

    #print("imax ", imax, " n ", n)
    for i in range(imax+1):
        if i > 0:
            #print("ideals ", ideals)
            ideals[i] = imul(ideals[i-1], ideals[i-1])
            if not ideals[i].is_reduced():
                ideals[i] = fast_reduce(ideals[i])
        if N - 2 ** (imax - i) >= 0:
            N -= 2 ** (imax - i)
        else:
            p[imax-i] = False
    sum = 0
    for i in range(len(p)):
        if p[i]:
           sum+= 2**i

    ideal_n = I(D)

    #ideal_n.print("ideal n")
    #ideals[0].print("0000")
    #fast_reduce(ideals[0]).print("1111")
    #fast_reduce(imul(ideal_n, ideals[i])).print("22222")
    for i in range(imax+1):
        #print("in range", p[i], len(p))
        if p[i]:
            ideal_n = fast_reduce(imul(ideal_n, ideals[i]))
            #ideal_n.print("ideal n")
    if n < 0:
        ideal0.b = - ideal0.b
    return ideal_n


def decryption(m):
    di = dict(zip([(ord(c) % 32) - 1 for c in string.ascii_letters], string.ascii_letters))
    base = 26
    step = 3
    messege = []
    for i in range(len(m)):
        reduce = m[i]
        for power in range(1, step +1):
            messege.append(di[reduce // base**(step - power)])
            reduce = reduce % base**(step - power)
    print(" messege ", messege, " len ", len(messege))




def imul(ideal1, ideal2):
    a1, b1 = ideal1.a, ideal1.b
    a2, b2 = ideal2.a, ideal2.b
    D = ideal1.D

    if((a1 == a2) and (b1 == b2)):
        d1, x1, y1 = a1, 1, 0
    else:
        d1, x1, y1 = egcd(a1, a2)
    a3 = a1 * a2
    d2, x2, y2 = egcd((b1 + b2)//2, d1)
    if d1 == d2:
        x2, y2 = 0, 1
    if(x1 == 0 and y1 !=0):
        x1, y1 = y1, x1
        b1, b2 = b2, b1
        a1, a2 = a2, a1
    a3 = a3//(d2 * d2)
    b3 = (y2 * x1 * (b2 - b1) + x2 * (D - b1 * b1)//(2 * a1)) % (2 * a2 // d2)
    b3 = b1 + b3 * (a1//d2) % (2 * a3)
    if(b3 > a3): b3 -= 2*a3
    if (b3 < -a3): b3 += 2 * a3
    # if abs(b3) > a3:
    #     b3 = b3 - 2 * a3 * nearest_int(b3, 2 * a3)
    return fast_reduce(ideal(a3, b3, D))
    #return (ideal(a3, b3, D))


def fast_reduce(idealA):
    a, b, D = idealA.a, idealA.b, idealA.D

    c = (b * b - D)//(4 * a)
   # print("b^2 ", b*b, " b^2-D ", b*b-D, " 4*a ", 4*a, " (b * b - D)/(4 * a) ", (b * b - D)/(4 * a))
   # print("c = ",c , " self c ", idealA.c)
    a, b, c = form_reduce((a, b, c))
    return ideal(a, b, D)

def form_reduce(form):
    form_normalized = form
    while not check_form_reduced(form_normalized):
        form_normalized = rho_form(form_normalized)
    return form_normalized

def rho_form(form):
    a, b, c = form
    s = nearest_int(b, 2 * c)
    return (c, -b + 2 * s * c, c * s * s - b * s + a)

def check_form_reduced(form):
    a, b, c = form
    if abs(b) > a or b == -a:
        return False

    if a > c or (a == c and b < 0):
        return False

    return True

def nearest_int(P, Q):
    q, R = P//Q, P % Q
    if 2 * R >= Q:
        return q + 1
    return q


# def primes(n):
#     """ Returns  a list of primes < n """
#     print("nnnnnnnnnnnnnn ", n, int(n**0.5)+1)
#     sieve = [True] * n
#     for i in range(3, int(n**0.5)+1, 2):
#         if sieve[i]:
#             print("int((n-i*i-1)/(2*i))+1 ", int((n-i*i-1)/(2*i))+1)
#             sieve[i*i::2*i] = [False]*(int((n-i*i-1)/(2*i))+1)
#     return [2] + [i for i in range(3, n, 2) if sieve[i]]


def construct_from_prime(idealA):

    print("Norm = ",idealA.norm())
    norm, bi, D = idealA.a, idealA.b, idealA.D
    #norm = fast_reduce(idealA).norm()
    prime_num = primes(norm + 1)
    iprime = construct_factor_base(idealA)
    omega_bar = [0] * (len(iprime) +1)
    F_temp = norm

    for i in range(len(iprime) -1, -1, -1):
        # print("i ", i , " iprime[i].norm() ",iprime[i].norm(), " F_temp ", F_temp)
        # print(" F_temp % (iprime[i].norm()) ", F_temp % (iprime[i].norm()))
        while F_temp % (iprime[i].norm()) == 0:
            F_temp = F_temp / iprime[i].norm()
            # print("here we are ", F_temp)
            omega_bar[i] += 1

    prip = []
    if F_temp == 1:
        omega_bar = omega_bar[:index_of_last_nonzero(lst = omega_bar)]
        # print("omega ", omega_bar)
    else:
        # print("norm cant be presented as prime ideal production. return (1,1) ideal"
        prip.append((I(D), 0))
        return I(D), prip

    for i in range(len(omega_bar)):
        if(omega_bar[i]):
            prip.append((iprime[i], omega_bar[i]))




    # print("prip ", prip)
    M  = 1
    for i in range(len(prip)):
        # print("prip[i][0].norm() ", prip[i][0].norm(), prip[i][1])
        M *= prip[i][0].norm()

    # print("M ", M)

    Mi = []
    Mir = []
    for i in range(len(prip)):
        Mi.append(M/prip[i][0].a)
        #print("i ", i, " Mi ", Mi[i])
    for i in range(len(prip)):
        # print("i ", i, Mi[i], prip[i][0].a, is_prime(prip[i][0].a))
        # print("Mi % a ", Mi[i]% prip[i][0].a )
        gcd_val, Mr, an = egcd(M/prip[i][0].a, prip[i][0].a)
        # print("i ", i, " Mi[i] ", Mi[i], " iprime[i].a ", prip[i][0].a, "gcd_val ", gcd_val, " Mr ", Mr, " an ", an)
        Mir.append(Mr)

    # for i in range(len(prip)):
    #     print("i ", i, " Mi*Mir ", Mi[i]*Mir[i]% prip[i][0].a, prip[i][0].a)
    # print("M",M, "mi ", Mi, "mir ", Mir)


    for i in range(2 ** len(prip)):
        ei = [1] * len(prip)
        binary = str(bin(i))[2:]
        count = 0
        while len(binary) < len(prip):
            binary = "0"+binary
            # print("binary ", binary)

        for echar in binary:
            # print("echar ", echar, " ",binary)
            if int(echar) == 0:
                ei[count] = -1
            else:
                ei[count] = 1
            count += 1
        b = 0
        for j in range(len(ei)):
            b += ei[j] * prip[j][0].b * Mi[j] * Mir[j]

        for j in range(len(ei)):
            if (bi - ei[j] * prip[j][0].b) % (2 * prip[j][0].a) == 0:
                # print("j ", j, " ei ", ei[j], " cond 11", (bi - ei[j] * prip[j][0].b) % (2 * prip[j][0].a) == 0)
                for k in range(len(prip)):
                    prip[k] = (prip[k][0], ei[k])
                ideal(norm, b, D).iprint("my i")
                fast_reduce(ideal(norm, b, D)).iprint("my fri")
                # for w in prip:
                #     print("prime ",w[0].a, " pow ", w[1])
                # return ideal(norm, b, D), prip
            #print("j ",j, " ei ", ei[j], " cond 11", (bi - ei[j] * prip[j][0].b) % (2 * prip[j][0].a) == 0)
            # print(" bi ", bi, " ei[j] * prip[j][0].b ", ei[j] * prip[j][0].b, "b - ei[j] * prip[j][0].b ", bi - ei[j] * prip[j][0].b)
            # print(" 2 * prip[j][0].a ", 2 * prip[j][0].a)

        #print("i ", i ," b ", b, " ei ", ei, cond)
        # ideal(norm, b, D).iprint("my i")
        # fast_reduce(ideal(norm, b, D)).iprint("my fri")



def index_of_last_nonzero(lst):
    for i, value in enumerate(reversed(lst)):
          if value != 0:
                 return len(lst) - i
    return len(lst)



def construct_factor_base(idealA):
    prime_numbers = primes(int(idealA.a) + 1)
    D = idealA.D
    iprime = []
    for i in range(len(prime_numbers)):
        for bp in range(-prime_numbers[i], prime_numbers[i] + 1):
            if (bp * bp) % (4*prime_numbers[i]) == D % (4*prime_numbers[i]):
                iprime.append(ideal(prime_numbers[i], bp, D))
                break
    return iprime


def prime_powers(start_ideal, factor_base):
    #iprime = construct_factor_base(start_ideal)
    iprime = factor_base
    D = start_ideal.D
    # while True:
    #     iprime = factorBase(start_ideal.D, 1)
    #     print("iprime ", iprime)
    #     if iprime != None:
    #         break
    # primes = []
    #
    # for i in range(len(iprime)):
    #     primes.append(iprime[i].a)




    #
    # prip = []
    # e = [0] * len(iprime)
    # delta_ideal = fast_reduce(start_ideal)
    # omega_bar = [0] * len(iprime)
    # for i in range(len(factor_base)):
    #     if delta_ideal.a == factor_base[i].a:
    #
    #         if fast_reduce(start_ideal).b >= 0:
    #             e[i] = 1
    #             omega_bar[i] = 1
    #         else:
    #             e[i] = -1
    #             omega_bar[i] = -1
    #         e = np.vstack((iprime, e))
    #         prip.append((iprime[i], omega_bar[i]))
    #         return delta_ideal, delta_ideal.a, prip, e, 1





    M = 5
    x_coord = 0
    completely_factors = False
    omega_bar = [0] * len(iprime)
    iter_count = 0
    while not completely_factors:
        print("M ", M, " iter_count ", iter_count)
        start_ideal.iprint("si")

        while True:
            e = np.random.randint(-1, 2, len(iprime))  # e = {-1, 0, 1}
            #print(" e ", e)
            e = np.vstack((iprime, e))
            #print("e ", e)
            delta_ideal = start_ideal

            for i in range(len(iprime)):
                #print("iprime ", iprime[i].a, iprime[i].b , " e ", e[1][i])
                delta_ideal = imul(delta_ideal, iexp(iprime[i], e[1][i]))  # d =  a*FB^e

            #delta_ideal.iprint("before")
            #delta_ideal = fast_reduce(delta_ideal)
            #delta_ideal.iprint("after")
            #print("power !!!!!!!!!!!!!", e)
            #delta_ideal.iprint("d = a*FB^e ")
            # print("delta_ideal norm", delta_ideal.norm(), " np.sqrt(abs(delta_ideal.D)/2)/M ", np.sqrt(abs(delta_ideal.D)/2)/M)
             #new line
            break
            #
            # if delta_ideal.norm() <= 0.5*np.sqrt(abs(delta_ideal.D)/2)/M or delta_ideal.norm() >= 5*np.sqrt(abs(delta_ideal.D)/2)/M:
            #     break




        #print("e ", e)
        for x in range(-M, M):
            omega_bar = [0] * len(iprime)
            F = delta_ideal.F(x)

            for i in range(len(factor_base)):
                if fast_reduce(delta_ideal).a == factor_base[i].a:
                    if fast_reduce(delta_ideal).b == factor_base[i].b:
                        b = factor_base[i].b
                    elif fast_reduce(delta_ideal).b == -factor_base[i].b:
                        b = -factor_base[i].b
                    F = factor_base[i].a
                    x = -(b + delta_ideal.b)/(2*delta_ideal.a)
                    # b = bA - 2 * x * v * aA - 2 * bA * v * y + 2 * cA * y * u
                    print("kskdskdksdkskdskdksdksdkskdksdksdkskdk")

            x_coord = x
            print(" F ", F)
            #print("1111x ", x, " F ", F, e)
            F_temp = int(F)
            # print("F_temp ", F_temp, len(iprime))

            for i in range(len(iprime)-1, -1, -1):
                # print("hereee ", i)
                while F_temp %  (iprime[i].norm()) == 0:
                    F_temp = F_temp/iprime[i].norm()
                    print(" F_temp ", F_temp)
                    omega_bar[i] += 1
            if F_temp == 1:
                completely_factors = True
                break

        iter_count += 1
        if iter_count % 100 == 0:
            M += 100



    F_test = 1
    for i in range(len(iprime) - 1, -1, -1):
        F_test *= iprime[i].norm()**omega_bar[i]
    # print("F_test ", F_test, F, " omega bar ", omega_bar, "F ", F, " e ", e)
    # delta_ideal.iprint("factor base")
    # fast_reduce(delta_ideal).iprint("reduce factor base")
    prip = []
    for i in range(len(omega_bar)):
        if(omega_bar[i]):
            prip.append((iprime[i], omega_bar[i]))

    return delta_ideal, F, prip, e, x_coord

def build_equivalent(idealA, f, x, y):
    aA = idealA.a
    bA = idealA.b
    cA = idealA.c
    D  = idealA.D
    a = f
    y = 1
    if y == 1:
        v = 1
        u = 0
    b = bA - 2*x*v*aA -2*bA*v*y +2*cA*y*u
    return ideal(a, b, D)

def construct_ideal_from_prime_powers(idealA, prip):
    if len(prip) == 0:
        return idealA, prip
    #print("in construct_ideal_from_prime_powers")
    a = idealA.a
    bi = idealA.b
    D = idealA.D
    #print("prip ", prip)
    M = 1
    for i in range(len(prip)):
        #print("prip[i][0].norm() ", prip[i][0].norm())
        M *= prip[i][0].norm()

    #print("M ", M)

    Mi = []
    Mir = []
    for i in range(len(prip)):
        Mi.append(M / prip[i][0].a)
        # print("i ", i, " Mi ", Mi[i])
    for i in range(len(prip)):
        #print("i ", i, Mi[i], prip[i][0].a, is_prime(prip[i][0].a))
        #print("Mi % a ", Mi[i] % prip[i][0].a)
        gcd_val, Mr, an = egcd(M / prip[i][0].a, prip[i][0].a)
        #print("i ", i, " Mi[i] ", Mi[i], " iprime[i].a ", prip[i][0].a, "gcd_val ", gcd_val, " Mr ", Mr, " an ", an)
        Mir.append(Mr)

    # for i in range(len(prip)):
    #     print("i ", i, " Mi*Mir ", Mi[i] * Mir[i] % prip[i][0].a, prip[i][0].a)
    # print("M", M, "mi ", Mi, "mir ", Mir)



    #print("len  prip ", prip)
    for i in range(2 ** len(prip)):
        ei = [1] * len(prip)
        binary = str(bin(i))[2:]

      #  print("lalalla binary ", binary)
        count = 0
        while len(binary) < len(prip):
            binary = "0"+binary
       #     print("binary ", binary)

        for echar in binary:
        #    print("echar ", echar, " ",binary)
            if int(echar) == 0:
                ei[count] = -1
            else:
                ei[count] = 1
            count += 1
        b = 0
        for j in range(len(ei)):
            if prip[j][0].b == prip[j][0].a:
         #       print("prip[j][0].b == prip[j][0].a ", prip[j][0].b == prip[j][0].a)
                prip[j][0].iprint(str(j))
                continue
            bp = (ei[j] * prip[j][0].b)% prip[j][0].a
            b += bp * Mi[j] * Mir[j]

        b = b%M
      #  print("la hula")
       # factorized_ideal = ideal(a, b, D)
       # print("al lula")
        # if not factorized_ideal.is_ideal():
        #     factorized_ideal.iprint("not ideal ideal(a, b, D)")
        #     # if int (factorized_ideal.c) != factorized_ideal.c:
        #     #     print("int(factorized_ideal.c)", int(factorized_ideal.c), " factorized_ideal.c ", factorized_ideal.c)
        #     #     if 4*int(factorized_ideal.c) != int(4*factorized_ideal.c):
        #     #         factorized_ideal.integer_mul(4)
        #     #         factorized_ideal.print("after mul ")
        #     #         ideal(factorized_ideal.a, factorized_ideal.b, factorized_ideal.D).print("my ideal after ")
        #     #     else:
        #     #         continue
        # else:
        #     factorized_ideal.iprint("norm ideal(a, b, D)")

        for j in range(len(ei)):
            p = prip[j][0].a
            bp = (ei[j] * prip[j][0].b)
            # if (b-bp)%p == 0:
            #     print("b? ",b, " p ", p," b%a ", b%p, " bp ", bp, "b%a = bp", True)
            # else:
            #     print("b? ", b, " p ", p, " b%a ", b % p, " bp ", bp, "b%a = bp", False)

        #print("binary ", binary, " b ", b, " bi ", bi, " M ",M)
        exitflag = True
        for j in range(len(ei)):
            #bp = (ei[j] * prip[j][0].b) % prip[j][0].a
            p = prip[j][0].a
            bp = (ei[j] * prip[j][0].b)# % prip[j][0].a
            # if (bi - bp) % (p) == 0:
            #     print("(bi - bp) % (a) == 0 ", j, (bi - bp) % (p) == 0, " bp ", bp)
            # else:
            #     print("(bi - bp) % (a) != 0 ", j, (bi - bp) % (p) != 0)

            if (bi - bp) % (2 * p) == 0:
                #print("j ", j, " bp ", bp, " a = ", p," exit  condition ", (bi - bp) % (2 * p) == 0)
                pass
            else:
                #print("j ", j, " bp ", bp, " a = ", p, " exit  condition ", (bi - bp) % (2 * p) == 0)
                exitflag = False
                continue

                # if not ideal(a, b, D).is_ideal():
                #     continue
        #print("i ", i, " e ", ei, " exit ", exitflag)
        if exitflag:
            for k in range(len(prip)):
                prip[k] = (prip[k][0], ei[k]*abs(prip[k][1]))

            result = I(D)

            for k in range(len(prip)):
                result = imul(result, iexp(prip[k][0], prip[k][1]))

            #result.iprint("result laldhdhdhd ")
            if not result.is_ideal():
                time.sleep(10)
         #   fast_reduce(result).iprint("fr result ")
         #   for w in prip:
         #      print("prime ",w[0].a, " pow ", w[1])
            return result, prip


def Lx(x, a, b):
    if x > math.e:
        return np.exp(np.log(np.log(x))**(1-a)*b*np.log(x)**(a))
    else:
        print("x < e return 0")
        return 0

def primeForm(D, p):
    #print("numberOfPrimeForms(D, p) ", numberOfPrimeForms(D, p))
    #if numberOfPrimeForms(D, p) > 0:
    if numberOfPrimeForms(D, p) > 0:
        #print("in primeForm")
        b = sqrtMod4P(D, p)
        #print("b ", b)
        if b is not None:
            return ideal(p, b, D)
    return None

def numberOfPrimeForms(D, p):
    return legendre_symbol(D, p) + 1

def sqrtMod4P(D, p):
    #print("legendre_symbol(D, p) ", legendre_symbol(D, p), legendre_symbol(p, D))
    if legendre_symbol(D, p) == 1:
        if p == 2:
            if D%2 == 0:
                return 2*(D/4 % 2)
            else:
                return 1
        else:
            r = sqrtModP(D, p)
            #print("r rrrrrrrrrrrrrr ", r, p, p-r)
            if r == None:
                return None
            if r%2 == D%2:
                return r
            else:
                return p-r
    else:
        #print("D is not square modulo p. returt 0")
        return None

def sqrtMod4PE(D, p, e):
    if D%p == 0:
        n = 1
        while 2*n+2 <= e and D%(p**(2*n+2)) == 0:
            n += 1
        if p ==2:
            number = modPower(2, -2*n, 4) # may be p must be prime
            number = number*D % 4
            if number!= 0 or number!=1:
                n = n -2
        else:
            D0 = D/p**(2*n)
            e0 = e -2*n
        if e0 == 0:
            return (p**n) * (D0%2)
        else:
            return (p ** n) * sqrtMod4P(D0, p**e0)
    else:
        b = sqrtMod4P(D, p)
        if b == None:
            return None
        f = 1
        while f < e:
            k = modPower(b, -1, p)
            k = k * (D - b*b)/(4*p**f)
            k = k %p
            b = (b + 2*k*p**f)%(2*p**(f+1))
            f += 1
        if b > p**e:
            return b - 2*p**e
        else:
            return b



def sqrtModP(r, p):
    m = p-1
    t = 0
    while m % 2 == 0:
        t += 1
        m /= 2
    #print("p ", p)
    #c = np.random.randint(1, p)
    for c in range(p+1):
        #print("ls ", legendre_symbol(c, p), legendre_symbol(p, c))
        if legendre_symbol(c, p) == 1:
            # if legendre_symbol(p, c) == 1:
            #print("sqrtModP fail")
            continue
        else:
            rm = modPower(r, m, p)
            cm = modPower(c, -m, p)
            e = 0
            i = 0
        while rm != 1:
            #print("here", rm, i, cm, m, c, r, p, t)
            #time.sleep(10)
            i += 1
            cm = cm * cm
            # print("here222",((rm**(2**(t-i-1)))%p))
            if ((rm ** (2 ** (t - i - 1))) % p) != 1:
                e += 2 ** i
                if rm == (rm * cm) % p:
                    break
                else:
                    rm = (rm * cm) % p
        if rm != 1:
            continue

        a = modPower(c, -e, p)
        a = (a*r)%p
        c_pow_tmp = 1
        for i in range(int(e/2)):
            c_pow_tmp *= c
            c_pow_tmp = c_pow_tmp % p

        a_pow_tmp = 1
        for i in range(int((m+1)/2)):
            a_pow_tmp *= a
            a_pow_tmp = a_pow_tmp % p
        return (c_pow_tmp*a_pow_tmp) % p
        #return (c**(e/2)*a**((m+1)/2))%p
    print("return None")
    return None


def modPower(num, pow, p):

    result = 1
    for i in range(int(abs(pow))):
        result *= num
        result = result % p
    #print(" (num**pow) % p ", result % p, num, pow, p)
    if pow >= 0:
        return result
    else:
        return egcd(result, p)[1]
        # return (num**pow) % p

        #return egcd(num**abs(pow), p)[1]


def kronecker(m, n):
    if n%2 == 0 and m%2 == 0:
        return 0
    if np.sign(n) == np.sign(m) and np.sign(n) == -1:
        j = -1
    else:
        j = 1
    n = abs(n)
    m = abs(m)
    if n == 0:
        if abs(m) == 1:
            return 0
        else:
            return 0
    while n%2 ==0:
        if modPower(m, 1, 8) == 3 or modPower(m, 1, 8) == 5:
            j = -j
        n = n/2
    m = modPower(m, 1, n)
    while m !=0:
        while m % 2 == 0:
            if modPower(n, 1, 8) == 3 or modPower(n, 1, 8) == 5:
                j = -j
            m = m / 2
        if modPower(m,1,4) ==3 and modPower(n,1,4) ==3:
            j = -j
        temp = m
        m = n
        n = temp
        m = modPower(m, 1, n)
    if n ==1:
        return j
    else:
        return 0

def primePowerForm(D, p , e):
    #b = sqrtMod4P(D, p) #hz
    b = sqrtMod4PE(D, p ,e)
    if b == None:
        return None
    c = (b*b - D)/(4*p**e)
    return ideal(p**e, b, c)


def factorBase(D, z):
    L = math.ceil(Lx(abs(D), 1/2, z))
    k = math.ceil(L)
    print("L ", L, " k ", k)
    F = []
    L = L + 1
    for p in primes(L):
        #print("kronecker(D, p) ", kronecker(D, p), legendre_symbol(D, p), p)
        #print("kronecker(p, D) ", kronecker(p, D), legendre_symbol(D, p), p)


        if kronecker(p, D) == 1:
        #if kronecker(D, p) == 1:
            print("p ", p)

            g = primeForm(D, p)
            if g == None:
                continue
            else:
                #print("b = ", g.b, " b delta p ", sqrtMod4P(D, p), " p ", (g.b-sqrtMod4P(D, p))%(2*p), " sigma ",(g.b+sqrtMod4P(D, p))%(2*p))
                F.append(g)
                #F.append(g.conjugation())
    return F

def randomRelation(FB, w):
    #print("lala")
    #w = [0] * len(fb)
    f = len(FB)
    if f == 0:
        return [None]
    D = FB[0].D
    # fb_with_conjugation = []
    # for prime_ideal in FB:
    #     fb_with_conjugation.append(prime_ideal)
    #     fb_with_conjugation.append(prime_ideal.conjugation())

    v = np.random.randint(0, abs(D), f)
    #print("00000v = ", v)
    ia = I(D)
    #print("v ", v)
    #print("w ", w)
    power = v + w

    #print("power ", power)
    for i in range(f):
        #print("i ",i)
        ia = imul(ia, iexp(FB[i], power[i]))
    #ia.iprint("1 ia ")
    ia = fast_reduce(ia)
    #bcfactor_time = datetime.now()
    cfactors, omega = completely_factors(ia, FB)
    #print("completely_factors time  ", datetime.now() - bcfactor_time)
    #print("cfactors ", cfactors)
    prip = []
    if cfactors:
        #ia.iprint("completely_factors ia")
        for i in range(len(omega)):
            if (omega[i]):
                prip.append((FB[i], omega[i]))
        cifpp_time = datetime.now()
        if len(prip) == 0:
            return [None]
        result, prip = construct_ideal_from_prime_powers(ia, prip)

      #  print("construct_ideal_from_prime_powers time  ", datetime.now() - cifpp_time)
        #ia.iprint("ia ")
        #fast_reduce(result).iprint("result ")
        z = power
     #   print("before z ", z)
        #print("len z ", len(z), " len(prip) ", len(prip), "omega ", omega)
        # for j in range(len(z)):
        #     for k in range(len(prip)):
        #         if prip[k][0].a == FB[j].a:
        #             z[j] -= prip[k][1]
        # print("after z ", z)
        for j in range(f):
            for k in range(len(prip)):
                #print("k ", k, "j", j, len(z), len(FB))
                #prip[k][0].iprint(str(j)+" "+str(k))
                if prip[k][0].a == FB[j].a:
                    if prip[k][0].b == FB[j].b:
                        z[j] -= prip[k][1]
      #  print("after z ", z)

        eideal = I(D)
        for k in range(len(z)):
       #     print("zk1111111111111 ", z[k])
            eideal = imul(eideal, iexp(FB[k], z[k]))
            #fast_reduce(eideal).iprint("eideal ")

        if all(zi == 0 for zi in z):
            return [None]
        return power
    else:
        #ia.iprint("ia ")
        return [None]

def fullRank(FB, z):

    f = len(FB)
    if f == 0:
        return [None]
    D = FB[0].D

   # print("f ", f)
    l = math.ceil(np.log(f)/p(D, z))
   # print("fllllllllllll ", l)
    B1 = math.ceil((f -1)*abs(D) + np.log(abs(D)))
    zf = []
    for i in range(f):
        print(" i ", i)
        e = [0] * f
        e[i] = 1
        e[i] = B1
       # print("e ", e)
        j = 0
        #zi = [None]*len(FB)
        while True:#j <= l :
            j += 1
            print("i", i, " j ", j)
           # print("e B1 ", B1)
            zi = randomRelation(FB, e)
            #print("in full rank ", zi)
            if None in zi:
                continue
                #return [None]
            else:

               # inds = sym.Matrix(np.matrix([z for z in zf])).T.rref()
               # print("inds before ", inds, " i ", i )

                zf.append(zi)
               # relation_matrix = np.matrix([z for z in zf])

                inds = sym.Matrix(np.matrix([z for z in zf])).T.rref()
                #print("inds after ", inds, " i ", i, " zf ", zf )
                print("len1 ", len(inds))
                if(len(inds) < i):
                    zf.pop()
                print("len2 ", len(inds))
             #   print("inds after pop", inds, " i ", i, " zf ", zf)
                # print("inds ", relation_matrix[list(inds[1])])
                # relation_matrix = relation_matrix[list(inds[1])]
                print("len(zf) ", len(zf), " i+1 ", i+1)
                if len(zf) == i+1:
                    print("break here")
                    break
        # if None not in zi:
        #     return [None] * len(FB)
        # else:
        #     zf.append(zi)
    relation_matrix = np.matrix([z for z in zf])
    #print("matrix ", relation_matrix)


    return relation_matrix


def p(D, z):
    return Lx(np.abs(D), 1/2, -1/(4*z) + g(D, z))

def g(D, z):
    return 1

def completely_factors(idealA, fb):
    omega_bar = [0] * len(fb)
    a = idealA.a
    #print("in cf a ", a, " len(fb) ", len(fb))
    for i in range(len(fb) - 1, -1, -1):
     #   print("i ",i, fb[i].norm())
        while a % (fb[i].norm()) == 0:
            a = a / fb[i].norm()
            omega_bar[i] += 1
      #      print("alalal ", a)
    if a == 1:
        a_test = 1
        for i in range(len(fb) - 1, 1, -1):
            a_test *= fb[i].norm() ** omega_bar[i]
       # print("a_test ", a_test, " a ", a, " omega bar ", omega_bar)
        return True, omega_bar
    else:
        return False, [0] * len(fb)

def egcd(a, b):
	x,y, u,v = 0,1, 1,0
	while a != 0:
		q, r = b//a, b % a
		#print("q = ",q ," r = ",r)
		m, n = x-u*q, y-v*q
		b,a, x,y, u,v = a,r, u,v, m,n
	gcd_val = b
	#print("x = ", x, " y = ", y)
	return gcd_val, x, y

def relationLattice(FB, z, mod):
    f = len(FB)
    if f == 0:
        return [None]
    D = FB[0].D
    N = 1
    i = 0

    B2 = (f+1)*abs(D) + np.log(abs(D))

    k = math.ceil(f*np.log2(B2))
    n = math.ceil(2**17 * np.log(k))
    l = math.ceil(np.log(n*k/p(D, z)))
    #print("additinal relation ", l," k ", k, " n ", n)
    z = []
    #return z
    Niter = k * n * l
    #Nmax = 2*int(np.sqrt(np.abs(D)))#k*n
    print("f ", f)
    if f < 10000000:
        Nmax = 550*f#abs(D)
    else:
        Nmax =  k*n

    while i < Niter and N <= Nmax:
        if i % 10 == 0:
            print("i = ",i, "Niter ", Niter, " Nmax ", Nmax, " N = ", N)
        #now = datetime.now()
        zi = randomRelation(FB, 0)
        #delta_time = datetime.now() - now

     #   print("N = ", N, " l ", l, " zi ", zi, " delta_time ", delta_time)
     #   print("additinal relation ", l, " k ", k, " n ", n)
        if None not in zi:
      #      print("zi all ")
            N += 1
            z.append(zi)
        i += 1
    if N <= Nmax:
        return [None]
    else:
        relation_lattice = np.matrix([zi%mod for zi in z])
        return relation_lattice



def power_representation(u):
    #print("u = ", u)
    N = u
    imax = 0
    while int(N / 2) != 0:
        imax += 1
        N = N / 2
    N = abs(u)
    p = [True] * (imax+1)
    for i in range(imax+1):
        if N - 2 ** (imax - i) >= 0:
            N -= 2 ** (imax - i)
        else:
            p[imax-i] = False
    sum = 0
    for i in range(len(p)):
        if p[i]:
           sum+= 2**i
    #print("sum ", sum, p)
    return p

# def reducePowerProduct(FB, u):
#     f = len(FB)
#     if f == 0:
#         return None
#     D = FB[0].D
#     reduced = ideal(1, 1, D)
#     for i in range(f):
#         print("i ", i)
#         if u[i] == 0:
#             q = ideal(1, 1, D)
#         else:
#             q = FB[i]
#         p = power_representation(u[i])
#         print("p ", p)
#         #for j in range(len(p)-1, -1, -1):
#         for j in range(len(p)):
#             print("j ", j)
#
#             if p[j]:
#                 if j == 0 :
#                     q2 = ideal(1, 1, D)
#                 else:
#                     q2 = fast_reduce(imul(q, q))
#             if p[j]:
#                 q = fast_reduce(imul(reduced, q2))
#         reduced = fast_reduce(imul(reduced, q))
#     return reduced



def find_li_vectors(dim, R):

    r = matrix_rank(R)
    index =  zeros( r ) #this will save the positions of the li columns in the matrix
    print("index ", index, " r ", r)
    counter = 0
    index[0] = 0 #without loss of generality we pick the first column as linearly independent
    j = 0 #therefore the second index is simply 0

    print("R.shape[0] ", range(r))
    for i in range(R.shape[0]): #loop over the columns
        print("R[:,i] ", R[:,i])
        if i != j: #if the two columns are not the same
            print("i ", i, " j ", j)
            inner_product = dot( R[:,i], R[:,j] ) #compute the scalar product
            print("inner_product ", inner_product)
            norm_i = norm(R[:,i]) #compute norms
            norm_j = norm(R[:,j])

            #inner product and the product of the norms are equal only if the two vectors are parallel
            #therefore we are looking for the ones which exhibit a difference which is bigger than a threshold
            if np.abs(inner_product - norm_j * norm_i) > 1e-4:
                counter += 1 #counter is incremented
                index[counter] = i #index is saved
                j = i #j is refreshed
            #do not forget to refresh j: otherwise you would compute only the vectors li with the first column!!
    print("here")
    R_independent = zeros((r, dim))
    print("R_independent ",R_independent, R[0,:])
    i = 0
    #now save everything in a new matrix
    while( i < r ):
        print("here2", R[index[i],:])
        R_independent[i,:] = R[index[i],:]
        i += 1

    return R_independent


def fliv():
    matrix = np.array(
        [
            [0, 1, 0, 0],
            [0, 0, 1, 0],
            [0, 1, 1, 0],
            [1, 0, 0, 1]
        ])

    print("det mat ",np.linalg.det(matrix))

    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[0]):
            if i != j:
                inner_product = np.inner(
                    matrix[:, i],
                    matrix[:, j]
                )
                norm_i = np.linalg.norm(matrix[:, i])
                norm_j = np.linalg.norm(matrix[:, j])

                print('I: ', matrix[:, i])
                print('J: ', matrix[:, j])
                print('Prod: ', inner_product)
                print('Norm i: ', norm_i)
                print('Norm j: ', norm_j)
                if np.abs(inner_product - norm_j * norm_i) < 1E-5:
                    print('Dependent')
                else:
                    print('Independent')


def from_pari_to_np(pari_matrix):
    pari = cypari2.Pari()
    col_size, row_size = pari.matsize(pari_matrix)
    list = []
    for i in range(col_size):
        temp_list = []
        for j in range(row_size):
            temp_list.append(float(pari_matrix[i][j]))
        list.append(temp_list)
    return np.matrix(list)

def from_np_to_pari(np_matrix):
    pari = cypari2.Pari()
    pari.allocatemem(10**9)
   # print("after alloc ")
    relation_list = []
    y_dim = len(np_matrix.tolist()[0])
    x_dim = len(np_matrix.tolist())

    #print("before cycle ")
    for i in range(len(np_matrix.tolist())):
        col = np_matrix.tolist()[i]
#        y_lattice_dimention = len(col)
        for j in range(len(col)):
            #print("i", i ," j ", j)
            el = col[j]
            relation_list.append(el)
    #print("relation_list ", relation_list)
    return pari.matrix(x_dim, y_dim, relation_list)

def representation_over_gens(idealA, FB, v, HNF):
    pari = cypari2.Pari()

    D = idealA.D
    f = len(FB)
    U, V, SNF = (pari.matsnf(from_np_to_pari(HNF), 1))
    SNF =  from_pari_to_np(SNF)
    #U = from_pari_to_np(U)
    m = []
    for i in range(f):
        m.append(float(SNF[i][i]))
    print("m ", m)
    print("SNF ", SNF)
    print("U ", U)
    print("V ", V)
    print("HNF ", HNF)
    print("v = ", v)

    result = I(D)
    for k in range(f):
        result = imul(result, iexp(FB[k], v[k]))

    fast_reduce(idealA).iprint("start ideal ")
    fast_reduce(result).iprint("prime_ideals_powers ")

    w = [0] * f
    a = [0] * f

    #U = [[]*f]*f
    i = k = f - 1
    print("1i ", i, " k ", k)
    while i >= 0:
        if HNF[i][i] == 1:
            for j in range(k):
                if i != j:
                    v[j] -= v[i]*HNF[j][i]
                else:
                    v[i] = 0
        i -= 1
    i = j = 0
    print("2i ", i, " k ", k)
    while i <= k:
        if HNF[i][i] != 1:
            w[j] = v[i]
            j += 1
        i += 1
        print("w ", w)

    for i in range(f):
        for j in range(f):
            a[i] += (w[j]*U[i][j]) % m[i]
    print("a11 ", a)

    for i in range(len(a)):
        a[i] = int(a[i])
    return a

def representation_over_FB(starti, fb):
    print("representation_over_FB")
    D = starti.D
    delta, F, prime_ideals_powers, e, x = prime_powers(starti, fb)
  #  delta.iprint("delta ")
  #  fast_reduce(delta).iprint("fr delta ")
  #  print("F ", F, " x ", x)
    cideal = build_equivalent(delta, F, x, 1)
   # cideal.iprint("cideal")
   # fast_reduce(cideal).iprint("fr cideal")

    #
    # if(fast_reduce(cideal).a != fast_reduce(delta).a or fast_reduce(cideal).b != fast_reduce(delta).b):
    #     print("alarm ")
    # else:
    #     print("all right ")


    #print("prime_ideals_powers ", prime_ideals_powers)
    if not prime_ideals_powers:
        print("prime_ideals_powers is empty")
        return [0]*len(fb)
    else:
        id, prime_ideals_powers = construct_ideal_from_prime_powers(cideal, prime_ideals_powers)
        # if (fast_reduce(cideal).a != fast_reduce(id).a or fast_reduce(cideal).b != fast_reduce(id).b):
        #     print("alarm ")
        # else:
        #     print("all right ")

    v = -e[1]
    for i in range(len(e[0])):
        e[0][i].iprint(str(i))
        for j in range(len(prime_ideals_powers)):
            if prime_ideals_powers[j][0].a == e[0][i].a:
                v[i] += prime_ideals_powers[j][1]


    print("v = ", v)
    result = I(D)
    for k in range(len(v)):
        e[0][k].iprint(str(k))
        #print("v[k ] ", v[k])
        result = imul(result, iexp(e[0][k], v[k]))
    fast_reduce(starti).iprint("starti")
    fast_reduce(result).iprint("result ")
    return v

def inverse_element(el, mod):
    if egcd(el, mod)[0] > 1:
        print("gcd ", egcd(el, mod)[0])
        return -1
    a = 1
    while a/el != int(a/el):
        a += mod
    return a/el

def classNumber(D):
    h = 0
    b = D % 2
    while b <= np.sqrt(np.abs(D)/3):
        A = (b*b - D)/4
        a = max(1, b)
        while a <= np.sqrt(A):
            c = A/a
            if A % a == 0:
                gcd = egcd(a, egcd(b, c)[0])[0]
                print("h ", h, " A ", A, " a ", a, " b ", b, " gcd ", gcd, " c ", c)
                if gcd == 1:
                    if b == 0 or a == b or a == c:
                        h += 1
                    else:
                        h += 2
            a += 1
        b += 2
    return h


def h_star(x, D):
    h = 1
    for p in primes(x):
        h *= 1/(1-legendre_symbol((D), p)/p)
       # print("hstar ",h, primes(x))
    return h

def B(x, D):
    pi_bar = 22/7
    B = 3*math.floor(np.sqrt(abs(D)))
    B *= h_star(x, D)
    B /= 4 *pi_bar
    return B

def determinant(matrix):
    pari = cypari2.Pari()
    hb = Decimal(hadamar_bound(matrix))
    print("hb ",hb, type(hb), Decimal((hb).sqrt()))
    c = math.floor(Decimal(hb).sqrt())
    prime_numbers = primes(c)
    m_len = len(matrix)
    prime_product = 1
    prime_iterator = 0
    while prime_product < c:
        prime_product *= prime_numbers[prime_iterator]
        prime_iterator += 1
    prime_numbers = prime_numbers[:prime_iterator]
    dets_mod_p = []
    for p in prime_numbers:
        tmp_matrix = pari.matrix(m_len, m_len, [0]*m_len*m_len)
        for i in range(m_len):
            for j in range(m_len):
                tmp_matrix[i][j] = matrix[i][j] % p
        det = pari.matdet(tmp_matrix) % p
        dets_mod_p.append(det)
    det_matrix = chrth(prime_numbers, dets_mod_p)
    print("det_matrix ", det_matrix)
    return det_matrix

def hadamar_bound(matrix):
    hb = 1
    print("len m ", len(matrix), matrix[10][10])
    for i in range(len(matrix)):
        row_sum2 = 0
        for j in range(len(matrix)):
            row_sum2 += matrix[i][j]*matrix[i][j]
        hb *= row_sum2
    hbt = 1
    for i in range(len(matrix)):
        col_sum2 = 0
        for j in range(len(matrix)):
            col_sum2 += matrix[j][i]*matrix[j][i]
        hbt *= col_sum2
    return int(min(hb, hbt))

def chrth(p, r):
    M  = 1
    for i in range(len(p)):
        M *= p[i]
    Mi = []
    Mir = []
    for i in range(len(p)):
        Mi.append(M / p[i])
    for i in range(len(p)):
        gcd_val, Mr, an = egcd(M/p[i], p[i])
        Mir.append(Mr)
    x = 0
    for i in range(len(p)):
        x += r[i]*Mi[i]*Mir[i]
    return x%M

def squarefree(n):
    print("n ", n)
    for i in range (2, math.ceil(sqrt(sqrt(n)))):
        if n%(i**2)==0:
            return False
    return True

def is_fundamental_discriminant(D):
    if D % 4 == 1:
        return True
    elif D % 4 == 0:
        return (not is_discriminant(D/4)) and squarefree((D/4)%4)
    else:
        return False

def is_discriminant(D):
    if D % 4 == 1 or D % 4 == 0:
        return True
    else:
        return False

def I(D):
    if D%4 == 1:
        return ideal(1, 1, D)
    else:
        return ideal(1, 0, D)

def calculate_determinant(full_rank_matrix, FB):
    pari = cypari2.Pari()
    detH = 1
    D = FB[0].D
    b_constatn = 100
    B_class_number = B(b_constatn, D)
    common_iterator = 0
    relation_iterator = 0

    magic_number = 10000000
    f = len(full_rank_matrix)
    while detH < B_class_number or detH > 2*B_class_number:
        #detH = 1
        tmp_det = detH
        zi = randomRelation(FB, 1)
        #print('zi ', zi)
        if None in zi:
            continue

        if relation_iterator >= f*f:#f*f*f*f: still magic
            relation_iterator = 1
            while True:
                rnd_i = np.random.randint(0, f)
                tmp_relation = copy.deepcopy(full_rank_matrix[rnd_i])
                full_rank_matrix[rnd_i] = zi
                new_det = abs(int(pari.matdet(from_np_to_pari(full_rank_matrix))))
                if new_det != 0:
                    detH = new_det
                    break
                else:
                    full_rank_matrix[rnd_i] = tmp_relation
            continue
        for i in range(len(full_rank_matrix)):

            tmp_relation = copy.deepcopy(full_rank_matrix[i])
            #print("i_string befor ", tmp_relation)
            full_rank_matrix[i] = zi
            #print("i_string full_rank_matrix ", full_rank_matrix, type(full_rank_matrix))
            det_lattice = abs(int(pari.matdet(from_np_to_pari(full_rank_matrix))))
            #det_lattice = abs(int(np.linalg.det(full_rank_matrix)))
            #print("after")

            full_rank_matrix[i] = tmp_relation
            if common_iterator == 0:
                detH = det_lattice
            detH = egcd(detH, det_lattice)[0]
            if detH == tmp_det:
                relation_iterator += 1
            #print("detH", detH, "det_lattice ", det_lattice)
            # print("npm I", npm.I)
            #
            # [H, U] = pari.mathnf(M_lattice, 4)
            #
            # print("H ", H, type(H))
            #
            # for i in range(len(H)):
            #     detH *= H[i, i]
            #print("det H", detH)
            common_iterator += 1
            if common_iterator > magic_number:
                return None
            #print("full_rank_matrix ", full_rank_matrix)
            #print("i_string ", tmp_relation)

    return detH


def fullRankRelationMatrix(FB, z):
    f = len(FB)
    if f == 0:
        return [None]
    D = FB[0].D
    l = math.ceil(np.log(f)/p(D, z))
    l = 1000*f
    zf = []
    counter = 0
    e = [0] * f
    #print("l ", l)
    while counter <= l:
        counter += 1
        if counter % 1 == 0:
            #print(" counter ", counter)
            pass

        zi = randomRelation(FB, e)

        #print("zi ", zi)
        if None in zi:
            continue
        else:
            zf.append(zi)
            if len(zf) % (2*f) == 0:
                #print("zf ", zf, len(zf), f)

                # tmp_mat = np.matrix([z for z in zf])
                # print("before inds")
                # inds = sym.Matrix(tmp_mat).T.rref()
                # print("after inds")
                tmp_det = 0
                zf_tmp = []
                while tmp_det == 0:
                    random.shuffle(zf)
                    zf_tmp = zf[:f]
                    tmp_mat = np.matrix([z for z in zf_tmp])
                    tmp_det =  abs(np.linalg.det(tmp_mat))
                    #print("lf af ", tmp_det)
                zf = zf_tmp
                #print("list(inds[1]) ", list(inds[1]))
                # tmp = []
                # for index in inds[1]:
                #     tmp.append(zf[index])
                # zf = tmp

                #print("zf ", zf, len(zf), f)
                if len(zf) == f:
                    break
    if counter >= l:
        return [None]
    relation_matrix = np.matrix([z for z in zf])
    #print("matrix ", relation_matrix)
    return relation_matrix


def random_search(relation_matrix, detMatrix):
    pari = cypari2.Pari()
    size = len(relation_matrix)
    while True:
        random_matrix = np.random.randint(-1, 2, size * size)
        # for i in range(len(random_matrix)):
        #         random_matrix[i] = random_matrix[i]*int(np.random.uniform(1, 1000))
        random_matrix = pari.matrix(size, size, random_matrix * detMatrix)
        print("random_matrix ", random_matrix)
        new_matrix = relation_matrix + random_matrix
        print("new_matrix ", new_matrix)
        det_new = abs(pari.matdet(new_matrix))
        det_rel = abs(pari.matdet(relation_matrix))
        print("relation_matrix ", relation_matrix, pari.matdet(relation_matrix))
        print("det_new ", det_new, "det_rel", det_rel)
        if det_new < det_rel:
            relation_matrix = new_matrix
        if det_new == detMatrix:
            break
    print("relation_matrix ", relation_matrix)
    return relation_matrix

def random_relation_matrix(file_csv, f, mod):
    matrix = file_csv.sample(n=f)
    matrix = matrix.drop(['Unnamed: 0'], axis =1)
    #print("matrix ", type(matrix))
    matrix = matrix.as_matrix()
    matrix = np.matrix([z % mod for z in matrix])
    # for i in range(f):
    #     for j in range(f):
    #         matrix[i][j] = matrix[i][j] % mod

    return matrix