from ideal import ideal, primeForm, p, random_relation_matrix, calculate_determinant, fullRankRelationMatrix, I, is_fundamental_discriminant, is_discriminant, squarefree, is_prime, kronecker, B, hadamar_bound, determinant, legendre_symbol, classNumber, inverse_element, representation_over_gens, representation_over_FB, starting_ideal, from_np_to_pari, convert_massage, from_pari_to_np, iexp, decryption, imul, fast_reduce,  egcd, construct_from_prime, construct_factor_base, prime_powers, build_equivalent, construct_ideal_from_prime_powers, fullRank, modPower, randomRelation, find_li_vectors#, reduce
import numpy as np
from datetime import datetime, timezone
import random
import math
import copy
import pandas as pd
D= -7*29*31*13*71*43*139*1601*16000000000000000000000000000000000000000000000000000000000000000001*1234567898765432183*19171323571
D= -119803744691274506954196478027666448305422107338661319930897882676033654173463236052636461146755017224699855500426777291347885152825740310470470188644025223079332669700563976732928049691985604910619
D= -143
#D = -974*4
D = -1513*4 #non factoring
#D = -1582*4
#D = -13727 #all right
#D = -18285*4 #non factoring
#D = -14033*4 #non factoring
#D = -14637*4 #non factoring
#D = -14722*4 #all right
#D = -19427 #all right
#D = -22965*4 #non factoring
#D = -19677*4 #non factoring
#D = -19919
import sympy as sym
#D= -55581980884703
#D = -23953*4
#D = -23910*4 # посчиталось все таки))
import ideal as il
import primesieve
import time
from numpy.linalg import inv
b_constant = 100
subexponentiality = 1
import cypari2
def main():
    pari = cypari2.Pari()

    discriminants_file = open('test_discriminants.txt', 'r')
    for item in discriminants_file:
        print("item ", int(item))
    return
    if not is_fundamental_discriminant(D):
        print("D is_fundamental_discriminant ", is_fundamental_discriminant(D))
        return


    m = convert_massage("happy thanksgiving")
    print("m = ",m)
    starti = starting_ideal(D)
    starti.iprint("start ideal = ")
    fast_reduce(starti).iprint("reduced start ideal = ")

#    return
    print("factor ", pari.factor(D))
    fb = il.factorBase(D, subexponentiality)
    for prime in fb:
        prime.iprint("prime ideal ")

    #return
    B_class_number = B(b_constant, D)
    print("B ", B_class_number)

   # return
   # return
    ##### begin here ##############################
   #  delta, F, prime_ideals_powers, e, x = prime_powers(starti, fb)
   #  delta.iprint("delta ")
   #  fast_reduce(delta).iprint("fr delta ")
   #  #return
   #  print("F ", F, " x ", x)
   #  cideal = build_equivalent(delta, F, x, 1)
   #  cideal.iprint("cideal")
   #  fast_reduce(cideal).iprint("fr cideal")
   #
   #  # rel = inverse_element(2, 10)
   #  # print("rel ", rel)
   # # return
   #  if(fast_reduce(cideal).a != fast_reduce(delta).a or fast_reduce(cideal).b != fast_reduce(delta).b):
   #      print("alarm ")
   #  else:
   #      print("all right ")
   #
   #  #return
   #  print("prime_ideals_powers ", prime_ideals_powers)
   #  if not prime_ideals_powers:
   #      print("prime_ideals_powers is empty")
   #  else:
   #      id, prime_ideals_powers = construct_ideal_from_prime_powers(cideal, prime_ideals_powers)
   #      if (fast_reduce(cideal).a != fast_reduce(id).a or fast_reduce(cideal).b != fast_reduce(id).b):
   #          print("alarm ")
   #      else:
   #          print("all right ")
   #
   #
   #  print("e = ", e[1])
   #  v = -e[1]
   #
   #  print("v here = ", v)
   #  for i in range(len(e[0])):
   #      e[0][i].iprint(str(i))
   #      for j in range(len(prime_ideals_powers)):
   #          if prime_ideals_powers[j][0].a == e[0][i].a:
   #
   #              v[i] += prime_ideals_powers[j][1]
   #
   #  print("prime_ideals_powers = ", prime_ideals_powers)
   #  print("e = ", -e[1])
   #  print("v = ", v)
   #  result = I(D)
   #  for k in range(len(v)):
   #      e[0][k].iprint(str(k))
   #      print("v[k ] ", v[k])
   #      result = imul(result, iexp(e[0][k], v[k]))
   #  fast_reduce(result).iprint("result ")
   #
   #
   #  gamma = delta.a* x + delta.b/2 + np.complex(np.sqrt(np.complex(delta.D)))
   #  print("gamma ", gamma, "a ",delta.a, " x ", x, " b ", delta.b)
   #  fast_reduce(starti).iprint("start ideal")
   #
   #  #return
    ##### end here ##############################
    discriminants_file = open('test_discriminants.txt', 'r')
    class_number = open('class_number.txt', 'w')
    max_counter = 100
    det = 1
    for item in discriminants_file:
        print("D = ", int(item))
        # continue
        fb = il.factorBase(int(item), subexponentiality)
        begin_time = datetime.now()

        counter = 0
        while counter < max_counter:


            cycle_count = 0
            while True:
                cycle_count += 1
                #print("count ", cycle_count)
                relation_matrix = fullRankRelationMatrix(fb, subexponentiality)
                if None not in relation_matrix:
                    break

            detMatrix = calculate_determinant(relation_matrix, fb)
            det = detMatrix
            counter += 1
            print("det ", detMatrix, int(item))
            #time.sleep(1)
            # print('relation_matrix', relation_matrix)
            # print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ok")
        end_time = datetime.now()
        run_time = (end_time - begin_time)/max_counter
        class_number.write(item + " " + str( det) + " " + str(run_time) + " " + str(max_counter) + " " + str(subexponentiality))
        #break



    return
    #rel = pd.read_csv('relation_list_test.csv')
    det  = tmp_det = 1
    #mat = random_relation_matrix(rel, 1000*len(e[0]), detMatrix)
    mat = il.relationLattice(e[0], 1, detMatrix)

    # x dim =
    # print(" len(np_matrix.tolist() ", len(mat.tolist()))
    # print("mat ", mat)
    # return
    #
    # relation_list = []
    # x_dim = y_dim = len(np_matrix.tolist())
    #
    # # print("before cycle ")
    # for i in range(len(np_matrix.tolist())):
    #     col = np_matrix.tolist()[i]
    #     #        y_lattice_dimention = len(col)
    #     for j in range(len(col)):
    #         # print("i", i ," j ", j)
    #         el = col[j]
    #         relation_list.append(el)
    # # print("relation_list ", relation_list)
    # return pari.matrix(x_dim, y_dim, relation_list)
    
    mat = from_np_to_pari(mat)

    print("mymatrix ", mat)

    pari.allocatemem(10 ** 9)
    [Ha, Ua] = pari.mathnf(mat, 1)
    print("H ", Ha)
    return
    #Ha = pari.matrix(2, 2, [5, 0, 0, 2])
    Sa = pari.matsnf(Ha, 1)
    print("llalalau ", Sa)


    #Ua = pari.matrix(2, 2, [2, -4, -5, 60])
    U_inv = pari.matsolve(Ua, 1)
    E = Ua*U_inv




    #Ha = pari.matrix(2, 2, [2, 0, 0, 5])
    #Ha = pari.matrix(2, 2, [10, 0, 0, 1])
    print("E ", E)

    print("uinv ", U_inv)
    print("H ", Ha)
    print("Ua ", Ua)
    #det = pari.matdet(Ha)
    #rank = pari.matrank(A)
    #print("determinant ", det)
    #print("rank ", rank)
    print("U ", Ua)
    return
    print("M*U = ", mat * Ua)
    Sa = pari.matsnf(Ha, 1)
    print("S = ", Sa[0]*Ha*Sa[1], len(Sa[0]))
    print("V ", Sa[1])
    print("u ", Sa)
    myi = iexp(ideal(2, 1, D), 8)
    myi = imul(myi, iexp(ideal(2, -1, D), 3))
    return




    while tmp_det != abs(detMatrix):
        mat = random_relation_matrix(rel, len(e[0]), detMatrix)
        mat = from_np_to_pari(mat)
        #mat = from_pari_to_np(mat)
        #print("mat ", mat, type(mat))
        #return
        old_det = abs(int(pari.matdet(mat)))
        tmp_det = 1
        while True:
            if old_det == tmp_det:
                break
            old_det = tmp_det
            for k in range(-1, 2, 2):
                detMatrix = k * abs(detMatrix)

                # for i in range(len(e[0])):
                #     gcd_tmp = int(mat[i][0])
                #     for j in range(len(e[0])):
                #         print("type ", type(gcd_tmp), type(int(mat[i][j])))
                #         gcd_tmp = egcd(gcd_tmp, int(mat[i][j]))[0]
                #     if gcd_tmp != 1:
                #         for j in range(len(e[0])):
                #             mat[i][j] = mat[i][j] / gcd_tmp
                #         break

                for i in range(len(mat)):
                    for j in range(len(e[0])):
                        #print("i ",i,"j", j, "k",k)
                        tmp_matrix = copy.deepcopy(mat)
                        rel_det = abs(int(pari.matdet(mat)))
                        tmp_matrix[i][j] += detMatrix

                        #print("mat ", mat)
                        #print("tmp_matrix ", tmp_matrix)
                        tmp_det = abs(int(pari.matdet(tmp_matrix)))
                        #print("tmp_det ", tmp_det, " rel_det ", rel_det)
                        time.sleep(1)

                        if tmp_det < rel_det:
                            mat = tmp_matrix

                        # Sa = pari.matsnf(mat, 1)
                        # print("Sa", Sa)
                        #print("rel_det ", rel_det)


        #det = abs(np.linalg.det(mat))
        print("det", rel_det)
    print("mat ", mat)
    return
















    #return
    rel_det = det = np.linalg.det(relation_matrix)
    #detpari = abs(int(pari.matdet(from_np_to_pari(relation_matrix))))
    print("det np ", det)
    tmp_det = det

    relation_matrix = from_np_to_pari(relation_matrix)
    for i in range(len(relation_matrix)):
        for j in range(len(relation_matrix)):
            relation_matrix[i][j] = relation_matrix[i][j] % detMatrix


    print("relation_matrix 222 ", relation_matrix, pari.matdet(relation_matrix))
    size = len(relation_matrix)
    #return
    # testmatrix = [18, 21, 34, 32, 8, 10, 3, 6, 9,
    # 1, 16, -15, 18, 5, 33, 1, -7, 24,
    # 0, 14, -26, 4, 5, 18, 3, -15, 35,
    # 21, -5, -18, 29, 8, 29, 3, -10, 30,
    # -2, -24, 21, 20, 32, 22, 29, -1, 6,
    # -7, 31, -12, 13, 7, 35, 4, -6, 34,
    # 6, 6, -7, 2, 2, 35, 35, 9, 27,
    # -25, 34, -22, 11, 33, 29, 25, -1, 26,
    # 2, 10, -12, 17, 26, 29, 35, -3, 33]
    # testmatrix = pari.matrix(size, size, testmatrix)
    # print("testmatrix ", testmatrix, pari.matdet(testmatrix))
    # relation_matrix = testmatrix
    f = len(fb)
    B2 = (f+1)*abs(D) + np.log(abs(D))

    k = math.ceil(f*np.log2(B2))
    n = math.ceil(2**17 * np.log(k))
    l = math.ceil(np.log(n*k/p(D, 1)))
    relation_counter  = 0
    relation_list = []
    while relation_counter < k*l*n:
        if relation_counter % 1000 == 0 and relation_counter != 0:
            pdobject = pd.DataFrame(relation_list)
            pdobject.to_csv('relation_list_test.csv', mode='a', header=False)
            relation_list = []
        zi = randomRelation(e[0], 1)
        print(relation_counter, k*l*n)
        if None in zi:
            continue
        else:
            relation_list.append(zi)
            relation_counter += 1

    return




    while tmp_det != detMatrix:

        zi = randomRelation(e[0], 1)
        f = len(e[0])
        print('zi ', zi)
        if None in zi:
            continue
        while True:
            rnd_i = np.random.randint(0, f)
            tmp_relation = copy.deepcopy(relation_matrix[rnd_i])
            relation_matrix[rnd_i] = zi
            new_det = abs(int(pari.matdet(relation_matrix)))
            if new_det != 0:
                detH = new_det
                break
            else:
                relation_matrix[rnd_i] = tmp_relation

        for k in range(-1, 2, 2):
            detMatrix = k* abs(detMatrix)

            for i in range(len(relation_matrix)):
                gcd_tmp = int(relation_matrix[i][0])
                for j in range(len(relation_matrix)):
                    print("type ", type(gcd_tmp), type(int(relation_matrix[i][j])))
                    gcd_tmp = egcd(gcd_tmp, int(relation_matrix[i][j]))[0]
                if gcd_tmp != 1:
                    for j in range(len(relation_matrix)):
                        relation_matrix[i][j] = relation_matrix[i][j] / gcd_tmp
                    break

            for i in range(len(relation_matrix)):
                for j in range(len(relation_matrix)):
                    tmp_matrix = copy.deepcopy(relation_matrix)
                    rel_det = abs(int(pari.matdet(relation_matrix)))
                    tmp_matrix[i][j] += detMatrix

                    print("relation_matrix ", relation_matrix)
                    print("tmp_matrix ", tmp_matrix)
                    tmp_det = abs(int(pari.matdet(tmp_matrix)))
                    print("tmp_det ", tmp_det, " rel_det ", rel_det)
                    # time.sleep(1)
                    if tmp_det < rel_det:
                        relation_matrix = tmp_matrix
                    else:
                        tmp_matrix = relation_matrix




        print("det ", tmp_det)
    print('relation_matrix', relation_matrix, np.linalg.det(relation_matrix))


    return
    for i in range(len(e[0])):
        e[0][i].iprint(str(i))

    print("fullRank 1 ", relation_matrix)
    det = np.linalg.det(relation_matrix)
    print("det ", det)
    return
    rmat = from_np_to_pari(relation_matrix)

    h0 = pari.matdet(rmat)
    #
    #
    #
    #
    #
    # A = pari.matrix(3, 3, [1, 2, 3, 4, -5, 6, 7, 8, 9])
    # determinant(A)
    # print("A ", A)
    # print("hb ", hadamar_bound(A))
    # det = pari.matdet(A)
    # print("det A ", det)


   # # return
   #  mylist = []
   #  for i in range(len(relation_matrix.tolist())):
   #      col = relation_matrix.tolist()[i]
   #      for j in range(len(col)):
   #          el = col[j]
   #          mylist.append(el)
   #
   #  M = pari.matrix(len(relation_matrix.tolist()), len(relation_matrix.tolist()), mylist)
   #
   #  [H, U] = pari.mathnf(M, 1)
   #  print("H ", H)
   #  print("U ", U)
   #  print("M = ", M)
   #  print("M*U = ", M*U)
   #
   #  S = pari.matsnf(H, 1)
   #  print("S = ", S[0]*H*S[1], len(S[0]))
    #return
    relation_pool = il.relationLattice(fb, 1)


    print("relation_pool ", relation_pool)
    #return
    relation_iterator = 0
    detH = 1
    while detH < B_class_number or detH > 2*B_class_number:
        #detH = 1
        relation_lattice = relation_pool[relation_iterator:]
        if len(relation_lattice) < len(fb):
            break
        # inds = sym.Matrix(relation_lattice).T.rref()
        # if (len(inds[1]) == 1):
        #     relation_lattice = min(relation_lattice)
        # print("relation_lattice ", relation_lattice, " len(inds) ", len(inds[1]))
        #
        # print("matrix ", relation_lattice)
        # print("inds ", inds[1], len(inds), type(inds[1]))
        # print("inds ", relation_lattice[list(inds[1])])

        # relation_lattice = relation_lattice[list(inds[1])]

        relation_lattice = relation_lattice[:len(fb)]

        M_lattice = from_np_to_pari(relation_lattice)
        print("M_lattice = ", M_lattice)

        # npm = from_pari_to_np(M_lattice)
        det_lattice = abs(int(pari.matdet(M_lattice)))
        if relation_iterator == 0:
            detH = det_lattice
        detH = egcd(detH, det_lattice )[0]
        # print("npm", npm)
        # print("npm I", npm.I)
        #
        # [H, U] = pari.mathnf(M_lattice, 4)
        #
        # print("H ", H, type(H))
        #
        # for i in range(len(H)):
        #     detH *= H[i, i]
        print("det H", detH)
        relation_iterator += 1

    return
    print("U ", U, type(U))
    print("M = ", M_lattice)
    print("M*U = ", M_lattice*U)
    #
    S = pari.matsnf(H, 1)
    print("S = ", S)
    Smith_nf =  from_pari_to_np(S[2])
    print("Smith_nf ", Smith_nf)
    print("S = ", S[0]*H*S[1], len(S[0]))
    U = from_pari_to_np(S[0])
    U_inv = np.matrix(U.I)
    print("uinv ", U_inv)
    return
    print("e 121212121212121 ", e)
    print("Uinv 121212121212121 ", U_inv)
    generators = []
    for i in range(len(relation_lattice)):
        generator = I(D)
        for j in range(len(relation_lattice)):
            print("U_inv[j] ", U_inv[j], i, j)
            print("U_inv[j][i] ", U_inv[j][i], i, j)
            print("U_inv[j][i][0] ", U_inv[j][i][0], i, j)
            generator = imul(generator, iexp(e[0][i], U_inv[j][i]))
        print("generators ")
        generator.iprint(str(i))
        generators.append(generator)

    for j in range(len(generators)):
        prime = I(D)
        for i in range(len(relation_matrix)):
            prime = imul(prime, iexp(generators[i], U[i][j]))
        print("primes ")
        prime.iprint(str(j))

    b_ideal = starting_ideal(D)
    print("aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa ")
    a = representation_over_gens(starti, fb, v, relation_lattice)
    print("amauaudef ", a , type(a))
    vb = representation_over_FB(b_ideal, fb)
    b = representation_over_gens(starti, fb, vb, relation_lattice)

    first_ideal = second_ideal = I(D)
    for i in range(len(generators)):
        first_ideal = imul(first_ideal, iexp(generators[i], a[i]))
    for i in range(len(generators)):
        second_ideal = imul(second_ideal, iexp(generators[i], b[i]))
    starti.iprint("a ideal")
    fast_reduce(starti).iprint("a ideal")
    first_ideal.iprint("first_ideal ")
    b_ideal.iprint("b ideal")
    fast_reduce(b_ideal).iprint("b ideal")
    second_ideal.iprint("second_ideal ")

    m = []
    for i in range(len(Smith_nf)):
        m.append(Smith_nf[i][i])
    x = 1


    for i in range(len(a)):
        print("a[i] ", a[i], type(a[i]))
        gcd = egcd(a[i], b[i])[0]
        if gcd > 1:
            a[i] = a[i] / gcd
            b[i] = b[i] / gcd
        if egcd(b[i], m[i])[0] > 1:
            x = -1
        else:
            b_inv = inverse_element(b[i], m[i])
            x = a[i] * b_inv % m[i]
            # x = m - x


    print("x  ", x, " a[i] ", a[i], " b[i] ", b[i])
    if x >= 0:
        iexp(second_ideal, x).iprint("b^x = a")
    else:
        print("no such solution")

    # A = pari.matrix(2, 2, [12, 2, 1, 1])
    # A = pari.matrix(1, 1, [10])
    # [Ha, Ua] = pari.mathnf(A, 1)
    # print("hahahllalalau ", Ha)
    # Sa = pari.matsnf(Ha, 1)
    # print("llalalau ", Sa)

    return

    A = pari.matrix(6, 2, [25 ,20 , 10, 80, 10, 0 , 0, 10, 20, 15, 80, 75])
    #
    #A = pari.matrix(2, 4, [137, 47, 96, 36, 2, 12, 30, 20])
    #A = pari.matrix(2, 9, [2, 10, 3, 13, 10, 0, 0, 10, 23, 13,  20, 15, 5, 10, 12, 2, 18, 8])
    A = pari.matrix(4, 2, [ 1, 9, 20, 10, 10, 20, 33, 13])
    A = pari.matrix(4, 2, [16, 6, 12, 2, 18, 8, 14, 4])
    A = pari.matrix(3, 2, [14, 4,  12, 2, 18, 8])   # A = pari.matrix(4, 2, [16, 6, 10, 0, 21, 31, 13, 3])
    A = pari.matrix(2, 2, [2, 4, -2, 6])
    A = pari.matrix(3, 2, [10, 0, 12, 22, 5, 15])
    A = pari.matrix(10, 2, [10, 0, 12, 22, 7, 17, 14, 4, 5, 15, 81, 21, 68, 8, 93, 73, 23, 13, 19, 9])
    A = pari.matrix(2, 2, [12, 2, 1, 1])
    det = pari.matdet(A)

    # A = pari.mattranspose(A)

    #A = pari.matrix(4, 2, [10, 0, 10, 5, 18, 8, 23, 3])
    # #A = pari.matrix(4, 4, [3, 3, 1, 4, 0, 1, 0, 0, 0, 0, 19, 16, 0,0,0,3])
    # A = pari.matrix(3, 4, [0,2,3,0,-5,3,-5,-5,4,3,-5,4])

    myi = ideal(1,1,D)
    # P = pari.mateigen(A)
    # P_inv = pari.matsolve(P, 1)
    # print("eigenvectors ", P, P_inv)



    iexp(ideal(2,-1,D), 5).iprint(9)
    print("A ", A)

    [Ha, Ua] = pari.mathnf(A, 1)
    print("H ", Ha)
    Ha = pari.matrix(2, 2, [5, 0, 0, 2])
    Sa = pari.matsnf(Ha, 1)
    print("llalalau ", Sa)


    Ua = pari.matrix(2, 2, [2, -4, -5, 60])
    U_inv = pari.matsolve(Ua, 1)
    E = Ua*U_inv




    #Ha = pari.matrix(2, 2, [2, 0, 0, 5])
    #Ha = pari.matrix(2, 2, [10, 0, 0, 1])
    print("E ", E)
    U_inv = U_inv*100
    print("uinv ", U_inv)
    print("H ", Ha)
    print("Ua ", Ua)
    #det = pari.matdet(Ha)
    #rank = pari.matrank(A)
    #print("determinant ", det)
    #print("rank ", rank)
    print("U ", Ua)
    print("M*U = ", A * Ua)
    Sa = pari.matsnf(Ha, 1)
    print("S = ", Sa[0]*Ha*Sa[1], len(Sa[0]))
    print("V ", Sa[1])
    print("u ", Sa)
    myi = iexp(ideal(2, 1, D), 8)
    myi = imul(myi, iexp(ideal(2, -1, D), 3))

    iexp(myi, 8).iprint("exp ")
    myi.iprint("myi")


    U_inv = pari.matsolve(Sa[0], 1)
    print("uinv ", U_inv)
    for i in range(len(relation_matrix)):
        generator = ideal(1,1,D)
        for j in range(len(relation_matrix)):
            generator = imul(generator, iexp(fb_with_conjugation[i], U_inv[j][i]))
        generator.iprint(str(i))
        generators.append(generator)

   # print("gen ", generators)

    for j in range(len(generators)):
        prime = ideal(1,1,D)
        for i in range(len(relation_matrix)):
            prime = imul(prime, iexp(generators[i], Sa[0][i][j]))
        prime.iprint(str(j))


    return
    # print("e len ", len(e[0]))
    # print("e ", e)
    # print("prime_ideals_powers ", prime_ideals_powers)
    # for i in range(len(e[0])):
    #     for j in range(len(prime_ideals_powers)):
    #         if e[0][i].a == prime_ideals_powers[j][0].a:
    #             e[1][i] = prime_ideals_powers[j][1] - e[1][i]
    #
    # print("e res ", e)
    #
    # equivalent_ideal = ideal(1,1,D)
    # for i in range(len(e[0])):
    #     equivalent_ideal = imul(equivalent_ideal, iexp(e[0][i], e[1][i]))  # d =  a*FB^e
    # equivalent_ideal.print("equivalent_ideal ")
    # fast_reduce(equivalent_ideal).print(" fr equivalent_ideal ")
   # v =
    # g, eg = construct_from_prime(delta)
    # g.print("g ")
    # print("eg ", eg)

    #return


#     completely_factors = False
#     xx = 0
#     n = 0
#     omega_bar = [0] * len(iprime)
#     F = 0
#     while not completely_factors:
#         e = np.random.randint(-1, 2, len(iprime))  # e = {-1, 0, 1}
#         #print("e ", e)
#         factor_base = starti
#
#         for i in range(len(iprime)):
#             factor_base = imul(factor_base, iexp(iprime[i], e[i]))  # d =  a*FB^e
#         factor_base.print("fbbbbbbb ")
#         print("e ", e)
#         M = 5
#
#         for x in range(-M, M):
#             omega_bar = [0] * len(iprime)
#             F = factor_base.F(x)
#             #print("x ", x, " F ", F, e)
#
#             F_temp = int(F)
#
#             for i in range(len(iprime) - 1, 1, -1):
#                 while F_temp %  (iprime[i].norm()) == 0:
#                     F_temp = F_temp/iprime[i].norm()
#                     omega_bar[i] += 1
#             if F_temp == 1:
#                 #print("x ", x, " yep")
#                 completely_factors = True
#                 n = F
#                 break
#
#
#     F_test = 1
#     for i in range(len(iprime) - 1, 1, -1):
#         F_test *= iprime[i].norm()**omega_bar[i]
#     print("F_test ", F_test, F, " omega bar ", omega_bar, "F ", F, " n ", n)
#     factor_base.print("factor base")
#     fast_reduce(factor_base).print("reduce factor base")
#
#     M_omega = 1
#     for i in range(len(iprime)):
#         if(omega_bar[i]):
#             M_omega *= iprime[i].norm()
#
#     M  = 1
#     for i in range(len(iprime)):
#         # print("i ", i)
#         M *= iprime[i].norm()
#
#     # print("M ", M)
#     # print("M ", (int(M)/int(2)))
#     # print("M ",(int(M) / int(2))*int(2), M == (int(M) / int(2))*int(2))
#     Mi = []
#     Mir = []
#     Mi_omega = []
#     Mir_omega = []
#     for i in range(len(iprime)):
#         # print("i ", i, " pr ", iprime[i].a)
#         Mi.append(M/iprime[i].a)
#         # print("i ", i, " Mi ", Mi[i])
#
#     iprimeM = []
#     for i in range(len(iprime)):
#         if (omega_bar[i]):
#             Mi_omega.append(M_omega/iprime[i].a)
#             iprimeM.append(iprime[i])
#
#
# #
# #
# #
#     for i in range(len(iprime)):
#         print("i ", i, Mi[i], iprime[i].a, is_prime(iprime[i].a))
#         print("Mi % a ", Mi[i]% iprime[i].a )
#         gcd_val, Mr, an = egcd(M/iprime[i].a, iprime[i].a)
#         print("i ", i, " Mi[i] ", Mi[i], " iprime[i].a ", iprime[i].a, "gcd_val ", gcd_val, " Mr ", Mr, " an ", an)
#         Mir.append(Mr)
#
#
#
#     for i in range(len(iprimeM)):
#         gcd_val, Mr, an = egcd(M_omega / iprimeM[i].a, iprimeM[i].a)
#         print("i ", i, " Mi[i] ", Mi_omega[i], " iprime[i].a ", iprimeM[i].a, "gcd_val ", gcd_val, " Mr ", Mr,                  " an ", an)
#         Mir_omega.append(Mr)
#
# #        print("gcd_val ", gcd_val)
#     for i in range(len(iprime)):
#         print("i " ,i," Mi*Mir ", Mi[i]*Mir[i]% iprime[i].a, iprime[i].a )
#     print("M",M, "mi ", Mi)
#
#
#     for i in range(len(iprimeM)):
#         print("i " ,i," Mi*Mir omega ", Mi_omega[i]*Mir_omega[i]% iprimeM[i].a, iprimeM[i].a )
#     print("M_omega ",M_omega, " mi _omega ", Mi_omega)
#
#     x = 0
#     for i in range(len(iprime)):
#         x += iprime[i].b*Mi[i]*Mir[i]
#     x = x%M
#     print("x ", x)
#
#     x_omega = 0
#     for i in range(len(iprimeM)):
#         x_omega += iprimeM[i].b*Mi_omega[i]*Mir_omega[i]
#     x_omega = x_omega%M_omega
#     print("x omega", x_omega)
#
#
# #
# #
# #     icycle = 0
#     omega = [0]*len(iprime)
#     vector = [] * len(omega_bar)
#     b_matrix = np.array([0]*(2**len(iprime)))
#     solution_count = 0
#     for i in range(len(b_matrix)):
#         ei = [1] * len(iprime)
#         count = 0
#         for echar in str(bin(i))[2:]:
#             if int(echar) == 0:
#                 ei[count] = -1
#             else:
#                 ei[count] = 1
#             count += 1
#         b = 0
#         for i in range(len(ei)):
#             b += ei[i]*iprime[i].b*Mi[i]*Mir[i]
#             #print("i ",i, " b ", b, " ei ", ei, " iprime[i] ", iprime[i].b, " Mi[i] ", Mi[i], " Mir[i] ", Mir[i])
#         find_solution = True
#         for i in range(1, len(ei)):
#             #print(i,"(b - ei[i]*iprime[i].b) % (2*iprime[i].a) ", (b - ei[i]*iprime[i].b) % (2*iprime[i].a))
#             if (b - ei[i]*iprime[i].b) % (2*iprime[i].a) != 0:
#                 find_solution = False
#                 #break
#
#         if find_solution:
#             solution_count += 1
#             print("we find a solution ", ei)
#             break
#     solution = iprime[0]
#
#     for i in range(1, len(ei)):
#         omega[i] = ei[i] * omega_bar[i]
#
#
#     vector = omega - e
#     for i in range(len(iprime) - 1, 1, -1):
#             temp = iexp(iprime[i], vector[i])
#             solution = imul(solution, temp)
#     solution.print("solution ")
#     fast_reduce(solution).print("solution ")
#     starti.print("start ideal = ")
#     fast_reduce(starti).print("start ideal reduce ")
#
#
#
#     print("vector ", vector, e, omega)
#
#     # for i in range(1, len(ei)):
#     #     omega[i] = ei[i]*omega_bar[i]
#     #     solution *= iexp(iprime[i], omega_bar[i])
#     # solution.print("solution")
#     # fast_reduce(starti).print("start ideal reduce ")
#
#     print("solution_count ", solution_count)




    # for e in range(-1, 1, 2):
    #     while icycle < len(iprime):
    #         e*iprime[icycle].b

    print("-----------------------------------------------------")
    a = [2, 3, 7]
    r = [1, 2, 6]

    M  = 1
    for i in range(len(a)):
        print("i ", i)
        M *= a[i]

    print("M ", M)
    Mi = []
    Mir = []
    for i in range(len(a)):
        Mi.append(M / a[i])
        print("i ", i, " Mi ", Mi[i], " pr ", a[i])


    for i in range(len(a)):
        gcd_val, Mr, an = egcd(M/a[i], a[i])
        print("i ", i, " Mi[i] ", Mi[i], " iprime[i].a ", a[i], "gcd_val ", gcd_val, " Mr ", Mr, " an ", an)
        Mir.append(Mr)

    for i in range(len(a)):
        print("Mi*Mir ", Mi[i]*Mir[i]% a[i], a[i])
    print("M",M, "mi ", Mi)

    x = 0
    for i in range(len(a)):
        x += r[i]*Mi[i]*Mir[i]
    x = x%M
    print("x ", x)

    for j in range(len(a)):
        if int(x % a[j]) == r[j]:
            print("x ", x, " a ", a[j], " x%a ", x % a[j], " r ", r[j], "x%a = r", True)
        else:
            print("x ", x, " a ", a[j], " x%a ", x % a[j], " r ", r[j], "x%a = r", False)
    # sti = ideal(2,1, -143)
    # idn = sti
    # for i in range(11):
    #     idn = imul(sti, idn)
    #     if not idn.is_reduced():
    #         idn = fast_reduce(idn)
    #     idn.print(str(i))
    print("5^36 ", 5**36%73)
    print("1% 73 ", 1%73, 72%73)

    print("betta 0 ",68**36 % 73)



    idealsX = [starti]*len(m)
    idealsY = [starti] * len(m)
    idealsXY = [starti] * len(m)
    idealsYX = [starti] * len(m)
    randomX = np.random.randint(100000, 10000000, len(m))
    randomY = np.random.randint(100000, 10000000, len(m))
    print("randomX ", randomX)
    print("randomY ", randomY)
    M = []
    for i in range(len(m)):
        idealsX[i] = iexp(starti, randomX[i])
        idealsY[i] = iexp(starti, randomY[i])
        idealsXY[i] = iexp(idealsX[i], randomY[i])
        idealsXY[i].print("XY "+str(i))
        M.append((m[i] + idealsXY[i].a, idealsY[i].a, idealsY[i].b))
    print("M ",M)

    for i in range(len(m)):
        idealsYX[i] = iexp(ideal(M[i][1], M[i][2], idealsYX[i].D), randomX[i])
        idealsYX[i].print("YX "+str(i))


    recovery_m = []
    for i in range(len(m)):
        recovery_m.append(M[i][0] - idealsYX[i].a)
    print("recovery m ", recovery_m)
    decryption(recovery_m)
    return


if __name__ == '__main__':
    error = main()
    if error is not None:
        print('Error occured:', error)


# def chinese_remainder(n, a):
#     sum = 0
#     prod = reduce(lambda a, b: a * b, n)
#
#     for n_i, a_i in zip(n, a):
#         p = prod / n_i
#         sum += a_i * mul_inv(p, n_i) * p
#     return sum % prod
#
#
# def mul_inv(a, b):
#     b0 = b
#     x0, x1 = 0, 1
#     if b == 1: return 1
#     while a > 1:
#         q = a / b
#         a, b = b, a % b
#         x0, x1 = x1 - q * x0, x0
#     if x1 < 0: x1 += b0
#     return x1
#
#
# if __name__ == '__main__':
#     n = [3, 5, 7]
#     a = [2, 3, 2]
#     print (chinese_remainder(n, a))




