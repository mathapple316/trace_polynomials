from sympy import *
import numpy as np
import itertools
from time import time
from sympy import degree as deg_of
import sys

################################################################################
# FILE LOGGING
################################################################################
''' print into the logfile'''


def fprint(*args):
    print(*args, file=fileinfo)
    return


################################################################################
# EPSILON RELATED
################################################################################

def eps_maker(dim):
    ''' Using generator, it makes all possible epsilons '''
    fprint("start eps_generator, dim: ", dim)
    print("-----------------------------------------------------------------------")
    ones = np.ones(dim, dtype=int)
    count = 0
    output = np.zeros(dim, dtype=int)
    count = count + 1
    yield output - ones

    while count < (3 ** dim):
        output[0] = output[0] + 1
        for j in range(dim - 1):
            if output[j] == 3:
                output[j] = 0
                output[j + 1] = output[j + 1] + 1
        count = count + 1
        yield output - ones
    fprint("finish eps_generator")
    if count != (3 ** dim):
        print("ERROR!!!")
    return


def is_alternating(eps):
    '''check whether given tuple is alternating'''
    count = 0
    first_val = 0
    first_val_flag = False
    for i in range(np.size(eps)):
        if eps[i] != 0:
            count = count + eps[i]
            if first_val_flag == False:
                first_val = eps[i]
                first_val_flag = True
            elif count != first_val and count != 0:
                return False
    if count == 0: 
        return True
    else: 
        return False

def degree(eps_input):
    '''compute the "degree" eps'''
   
    eps = eps_input.copy()
    size = np.size(eps)
    sgn = 1
    count = 0
    subcount = 0
    # flag denotes eps[i] is in a subinterval
    flag = True
    # sgn_flag denotes next zero should change sgn.
    sgn_flag = False

    # Roll left until fir entry becomes 0
    eps = roll_until(eps, 0)

    # Cut out
    eps = cut_out(eps)

    # Roll left until first entry becomes 3
    eps = roll_until(eps, 3)

    for i in range(size):
        # start of subinterval
        fprint("eps[i] : ", eps[i], ", ", flag, sgn_flag)
        if eps[i] == 0 and not flag:
            sgn = 1
            subcount = 1
            flag = True
            sgn_flag = False

        # end of subinterval, or inside [3,3,...,3]
        elif eps[i] == 3:
            count = count + abs(subcount)
            fprint("count : ", count)
            subcount = 0
            flag = False
        # zero, inside a subinterval
        # since sgn can change after 0, reset sgn_flag
        elif eps[i] == 0 and flag == True:
            if sgn_flag:
                sgn = -sgn
                fprint("change sgn")
                sgn_flag = False
            subcount = subcount + sgn
        # nonzero, inside a subinterval
        elif eps[i] != 0 and flag:
            sgn_flag = True

        # last entry
        if i == size - 1:
            count = count + abs(subcount)
            fprint("count : ", count)

    return (int)(1 / 2 * count)


################################################################################
# MANUAL EXPANSION
################################################################################

def tr_manual(m):
    # initiation

    if np.size(m) % 2 != 0:
        dim = np.size(m) + 1
    else:
        dim = np.size(m)

    r = np.size(m) / 2

    if dim == 2:
        if np.all(m == 0):
            fprint("tr_for_basis, ", m, "return 2")
            return 2
        elif m[0] == 1 and m[1] == 0:
            fprint("tr_for_basis, ", m, "return x")
            return x
        elif m[0] == 0 and m[1] == 1:
            fprint("tr_for_basis, ", m, "return y")
            return y
        else:
            fprint("tr_for_basis, ", m, "return z")
            return z

    result = 0
    for idx in idx_generator(dim, 2):
        fprint("in for loop, idx:", idx, "dim: ", dim)
        summand = ((-1) ** sum(idx)) * (vector_trace(m, tr_for_idx(dim, idx)))
        result += summand
        fprint("result : ", result)
    return result


def tr_for_idx(dim, index):
    if dim % 2 != 0:
        fprint("ERROR! odd dimension")
        return

    if dim == 2:
        if np.all(index == 0):
            fprint("tr_for_idx, ", index, "return 2")
            return 2
        elif index == [1, 0]:
            fprint("tr_for_idx, ", index, "return x")
            return x
        elif index == [0, 1]:
            fprint("tr_for_idx, ", index, "return y")
            return y
        else:
            fprint("tr_for_idx, ", index, "return z")
            return z

    if np.all(index == 1):
        return beta(dim / 2, z)
    else:
        for i in range(0, dim):
            if index[i] == 0:
                reduced_idx = reduce_idx(index, dim, i)
                return tr_manual(reduced_idx)
            else:
                continue


def reduce_idx(index, dim, i):
    if dim == 2:
        fprint("ERROR!!!!")
        return

    reduced = np.zeros(dim - 2, dtype=int)

    '''iê°€ ì²˜ìŒ?— ?žˆ?Š” ê²½ìš°'''
    if i == 0:
        for j in range(0, dim - 2):
            if j == (dim - 3):
                reduced[j] = index[dim - 1] + index[1]
            else:
                reduced[j] = index[j + 2]

    elif i == (dim - 1):
        for j in range(0, dim - 2):
            if j == 0:
                reduced[j] = index[0] + index[dim - 2]
            else:
                reduced[j] = index[j]

    else:
        for j in range(0, dim - 2):
            if j < (i - 1):
                reduced[j] = index[j]
            elif j == (i - 1):
                reduced[j] = index[i - 1] + index[i + 1]
            else:
                reduced[j] = index[j + 2]
    return reduced


def idx_generator(dim, max_val):
    ''' Using generator, it makes all possible epsilons '''
    fprint("start index generator, dim: ", dim)
    count = 0
    output = np.zeros(dim, dtype=int)
    count = count + 1
    fprint("1 idx:", output, ", count:", count, "/", max_val ** dim)
    yield output

    while count < (max_val ** dim):
        output[0] = output[0] + 1
        for j in range(dim - 1):
            if output[j] == max_val:
                output[j] = 0
                output[j + 1] = output[j + 1] + 1
        count = count + 1
        fprint("2 idx:", output, ", count:", count)
        yield output

    fprint("finish idx_generator")
    if count != (max_val ** dim):
        fprint("ERROR!!!")

    return


def vector_trace(vector, word):
    return vector * tr_manual(word)


def vector_fullpoly(vector, full_poly):
    deg_z = deg_of(full_poly, z)
    result = 0

    for k in range(0, deg_z + 1):

        if coeff_of_kpower(full_poly, k) == 0:
            continue

        result += vector_coeffpoly(vector, coeff_of_kpower(full_poly, k)) * (z ** k)

    return expand(result)


# z^k?˜ coeffê°€ ?™”?‹¤ê³? ê°€? •
def vector_coeffpoly(vector, coeff_poly):
    sum = 0

    if deg_of(coeff_poly, x) == 0 and deg_of(coeff_poly, y) == 0:
        return coeff_poly * large_b(vector, x, y)
    else:
        summand_array = str(coeff_poly).split('+')
        num_of_summand = len(summand_array)

    for i in range(0, num_of_summand):
        sum += vector_coeffmono(vector, summand_array[i], np.zeros(len(vector), dtype=int))
    return sum


def vector_coeffmono(vector, coeff_mono, flag):
    x_deg = deg_of(coeff_mono, x)
    y_deg = deg_of(coeff_mono, y)

    if len(vector) / 2 < x_deg or len(vector) / 2 < y_deg:
        print("error!!!!!!")
        return

    if x_deg != 0:
        for i in range(0, x_deg):
            flag[2 * i] = 1

    if y_deg != 0:
        for i in range(0, y_deg):
            flag[2 * i + 1] = 1

    return large_b_with_flags(vector, flag)


def large_b_with_flags(vector, flag):
    if np.all(flag == 0):
        return large_b(vector, x, y)

    else:
        for i in range(0, len(flag)):
            if flag[i] == 1:
                flag[i] = 0
                vector1 = np.copy(vector)
                vector2 = np.copy(vector)
                vector1[i] = vector[i] + 1
                vector2[i] = vector[i] - 1
                return large_b_with_flags(vector1, flag) - large_b_with_flags(vector2, flag)
        print("Never reach!!!!!!")
        return


def coeff_of_kpower(expr, power):
    if power == 0:
        return expand(expr).as_coeff_add(z)[0]
    else:
        return expand(expr).coeff(z ** power)


################################################################################
# OPERATIONS ON EPSILON
################################################################################


def roll_until(eps, num):
    ''' roll(-1) [eps] until first entry becomes [num] '''
    size_eps = np.size(eps)

    if np.all(eps != num) == True:
        fprint("no entry is ", num, " stop rolling. ")
        return eps

    while (eps[0] != num):
        eps = np.roll(eps, -1)
    fprint('rolling: ', eps)
    return eps


def cut_out(eps_input):
    '''compute the maximum of the length of consecutive non-zero subintervals'''
    '''change the maximal non-zeros subintervals into [3,3,...3]'''
    '''3 will be used as a delimeter between subintervals '''
    # Find the maximum even-length of non-zero subintervals
    eps = eps_input
    size = np.size(eps)
    even_leng = size
    count = 0

    while (even_leng > 0):
        count = 0
        for i in range(np.size(eps)):
            if eps[i] == 3:
                count = 0
            else:
                count = count + abs(eps[i])
                if eps[i] == 0:
                    if count == even_leng:
                        for j in range(count):
                            eps[i - j - 1] = 3
                            fprint("1 middle cut out : ", eps, "even leng: ",
                                   even_leng, "i: ", i)
                            count = 0
                    else:
                        count = 0
                elif i == size - 1:
                    if count == even_leng:
                        for j in range(count):
                            eps[i - j] = 3
                            fprint("2 middle cut out : ", eps, "even leng: ",
                                   even_leng, "i: ", i)
                            count = 0
        even_leng = even_leng - 2
    fprint("cut out : ", eps)
    return eps


################################################################################
# Deprecated
################################################################################

def degree_of_subinterval(sub_eps):
    ''' compute "degree" of eps whose all maximal-even-length-subintervals
    were cut-out'''
    size = np.size(sub_eps)
    sign = 1
    count = 0
    flag = False
    if (sub_eps[0] != 0) or (sub_eps[size - 1] != 0):
        fprint("epsilon format error!!!!: ", sub_eps)
        return 0
    if np.size(sub_eps) == 1:
        return 0
    elif np.size(sub_eps) == 2:
        return 1
    for i in range(size):
        if sub_eps[i] == 0:
            count = count + sign
            flag = false
            fprint("i: ", i, "count: ", count)
        elif sub_eps[i] != 0 and not flag:
            sign = -sign
            flag = true
            fprint("i: ", i, "change sign")

    fprint(" degree : ", abs(int(count / 2)))
    return abs(int(count / 2))


def nchoosek(n, k):
    numerator = 1
    denominator = 1
    k = min(n - k, k)
    for i in range(1, k + 1):
        denominator *= i
        numerator *= n + 1 - i
    return int(numerator / denominator)


def beta_explicit(num, symbol):
    if num == 0:
        return 0
    elif num < 0:
        return (-1) * beta(-num, symbol)
    count = (num - 1) / 2
    total = 0
    for i in range(0, count + 1):
        total += ((-1) ** i) * nchoosek(num - 1 - i, i) * (symbol ** (num - 1 - 2 * i))
    return total


################################################################################
# RECURRENCE RELATIONS
################################################################################
def beta(num, symbol):
    if num < 0:
        return expand(-beta(-num, symbol))
    if num == 0:
        return 0
    elif num == 1:
        return 1
    else:
        ret = (symbol * beta(num - 1, symbol)) - beta(num - 2, symbol)
        return expand(ret)


def tau(num, symbol):
    if num < 0:
        return tau(-num, symbol)
    if num == 0:
        return 2
    elif num == 1:
        return symbol
    else:
        ret = (symbol * tau(num - 1, symbol) - tau(num - 2, symbol))
        return ret


def large_b(m, symbol1, symbol2):
    if np.any(m == 0):
        return 0
    fprint("large_b, input :", m)
    expr = 1
    multiplicand = 1
    for i in range(np.size(m)):
        if i % 2 == 0:
            multiplicand = beta(m[i], symbol1)
            fprint("i:", i, ", beta :", multiplicand)
            expr = expr * multiplicand
        else:
            multiplicand = beta(m[i], symbol2)
            fprint("i:", i, ", beta :", multiplicand)
            expr = expr * multiplicand
    return expand(expr)


################################################################################
# MAIN LOGIC
################################################################################

def inputAndInit():
    # input vector m and set dim
    print("\n type the value of m ")
    line = input(" ").strip().split(' ')
    print("\n m : ", line)

    if len(line) % 2 != 0:
        dim = len(line) + 1
    else:
        dim = len(line)

    r = (int)(dim / 2)

    global filename
    filename = "trace"

    m = np.zeros(dim, dtype=int)
    for i in range(len(line)):
        m[i] = (int)(line[i])
        filename = filename + "_" + str(line[i])

    filename = filename + ".txt"

    global fileinfo
    fileinfo = open(filename, 'w', -1, "utf-8")

    fprint("m : ", line)

    global z, x, y
    z, x, y = symbols('z,x,y')
    init_printing(order='rev-lex')

    return m


def tr(m):
    # initiation

    if np.size(m) % 2 != 0:
        dim = np.size(m) + 1
    else:
        dim = np.size(m)

    r = (int)(np.size(m) / 2)

    # Compute expr1
    expr1 = expand((1 / 2) * tau(r, z) * large_b(m, x, y))

    # Compute expr2
    gen = eps_maker(dim)
    alt_count = 0
    expr2 = 0
    summand = 0
    deg = 0
    largeb = 0
    for mu in gen:
        if is_alternating(mu) == True:
            alt_count = alt_count + 1
            fprint(alt_count, ". mu:", mu)
            deg = degree(mu)
            fprint(" degree of mu : ", deg, "\n")
            largeb = large_b(m + mu, x, y)
            if largeb == 0:
                summand = 0
                fprint("-----------------------------------------------------------------------")
                continue
            else:
                summand = expand(((-1) ** (r - deg)) * tau(deg, z) * largeb)
            fprint(" large_b : ", largeb)
            fprint(" sgn : ", ((-1) ** (r - deg)))
            fprint(" tau(deg) : ", tau(deg, z))
            fprint(" summand : ", summand)
            expr2 = expr2 + summand
            fprint("-----------------------------------------------------------------------")

    if alt_count != (2 ** dim) - 1:
        print("ERROR!!!")
    fprint(" expr1 : ", expr1)
    fprint(" expr2 : ", expr2)
    fprint("-----------------------------------------------------------------------")
    expr = expand(expr1 + (1 / 2) * expr2)
    fprint(" RESULT is : ")
    fprint(expr)
    # print("expr :", expr)
    # pretty_print(expr, order='rev-lex')
    return expr


def check(m, expr):
    if np.size(m) % 2 != 0:
        dim = np.size(m) + 1
    else:
        dim = np.size(m)

    A = Matrix([[2, 1], [0, 1 / 2]])
    B = Matrix([[3, 0], [1, 1 / 3]])
    C = Matrix([[1, 0], [0, 1]])

    for i in range(dim):
        if i % 2 == 0:
            C = C * A ** (m[i])
        else:
            C = C * B ** (m[i])

    fprint("-----------------------------------------------------------------------")
    fprint(" Calculation Example ")
    fprint(" A : ", A)
    fprint(" B : ", B)
    fprint(" trace(C) : ", trace(C))
    fprint(" expr.sub : ", expr.subs({x: trace(A), y: trace(B), z: trace(A * B)}))

    print("")
    print(" RESULT is : \n ")
    pretty_print(expr, order='rev-lex')
    print("\n Read " + filename + " for detailed calculation. ")


################################################################################

if __name__ == "__main__":
    sys.setrecursionlimit(10000)
    m = inputAndInit()
    t0 = time()
    expr = tr_manual(m)
    print("calculation time : ", time() - t0)
# check(m,expr)
