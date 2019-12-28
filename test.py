from sympy import * 
import numpy as np
import sys
from time import sleep

################################################################################
# DEBUGGER
################################################################################

def dprint(*args_tuple):
    ''' detailed print for debugging'''
    return

################################################################################
# EPSILON RELATED
################################################################################

def eps_maker(dim):
    ''' Using generator, it makes all possible epsilons '''
    dprint("start eps_generator, dim: ",dim)
    ones = np.ones(dim, dtype=int)
    count = 0
    output = np.zeros(dim, dtype=int)
    count = count+1 
    yield output - ones

    while(count < (3**dim)):
        output[0] = output[0]+1
        for j in range(dim - 1): 
            if output[j] == 3:
                output[j] = 0
                output[j+1] = output[j+1] + 1
        count = count + 1
        yield output - ones 
    dprint("finish eps_generator")
    if count != (3**dim):
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
    eps = roll_until(eps,0)

    # Cut out
    eps = cut_outt(eps)
    
    # Roll left until first entry becomes 3
    eps = roll_until(eps, 3) 

    for i in range(size):
        # start of subinterval
        dprint("eps[i] : ", eps[i],", ", flag, sgn_flag)
        if eps[i] == 0 and flag == False:
            sgn = 1
            subcount = 1
            flag = True
            sgn_flag = False
                
        # end of subinterval, or inside [3,3,...,3]
        elif eps[i] == 3 :
            count = count + abs(subcount)
            dprint("count : ", count)
            subcount = 0
            flag = False
        # zero, inside a subinterval
        # since sgn can change after 0, reset sgn_flag
        elif eps[i] == 0 and flag == True:
            if sgn_flag == True:
                sgn = -sgn
                dprint("change sgn")
                sgn_flag = False
            subcount = subcount + sgn
        # nonzero, inside a subinterval
        elif eps[i] != 0 and flag == True:
            sgn_flag = True

        # last entry
        if i == size - 1:
                count = count + abs(subcount)
                dprint("count : ", count)

    return (int)(1/2 * count)

################################################################################
# OPERATIONS ON EPSILON
################################################################################


def roll_until(eps, num):
    ''' roll(-1) [eps] until first entry becomes [num] '''
    size_eps = np.size(eps) 

    if np.all(eps != num) == True:
        dprint("no entry is ", num," stop rolling. ")
        return eps

    while (eps[0] != num):
        eps = np.roll(eps,-1) 
    dprint('rolling: ',eps)
    return eps

def cut_out(eps_input):
    '''compute the maximum of the length of consecutive non-zero subintervals'''
    '''change the maximal non-zeros subintervals into [3,3,...3]'''
    '''3 will be used as a delimeter between subintervals '''
    # Find the maximum even-length of non-zero subintervals
    eps = eps_input
    max_leng = 0
    count = 0
    size = np.size(eps)
    for i in range(size):
        if eps[i] != 0:
            count = count + 1
        if eps[i] == 0 or i == size-1:
            if (count%2 == 0 and count > max_leng):
                max_leng = count
            count = 0
    for i in range(np.size(eps)):
        count = count + abs(eps[i])
        # 0 entries
        if eps[i] == 0: 
            if (count%2 == 0 and count == max_leng):
                for j in range(1,count+1):
                    eps[i-j] = 3
            count = 0
        elif i == size-1:
            if (count > 0 and count%2 == 0 and count == max_leng):
                for j in range(count):
                    eps[i-j] = 3
            count = 0
            dprint("cut out : ", eps)
    return eps

def cut_outt(eps_input):
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
                            eps[i-j-1] = 3
                            dprint("1 middle cut out : ", eps, "even leng: ", even_leng, "i: ", i)
                            count = 0
                    else:
                        count = 0
                elif i == size - 1:
                    if count == even_leng:
                        for j in range(count):
                            eps[i-j] = 3
                            dprint("2 middle cut out : ", eps, "even leng: ", even_leng, "i: ", i)
                            count = 0
        even_leng = even_leng - 2
    dprint("cut out : ", eps)
    return eps
################################################################################
# deprecated
################################################################################

def degree_of_subinterval(sub_eps):
    ''' compute "degree" of eps whose all maximal-even-length-subintervals 
    were cut-out'''
    size = np.size(sub_eps)
    sign = 1
    count = 0
    flag = False
    if (sub_eps[0] != 0) or (sub_eps[size - 1] != 0):
        dprint("epsilon format error!!!!: ", sub_eps)
        return 0
    if np.size(sub_eps) == 1:
        return 0
    elif np.size(sub_eps) == 2:
        return 1
    for i in range(size):
        if sub_eps[i] == 0:
            count = count + sign
            flag = false
            dprint("i: ",i,"count: ",count)
        elif (sub_eps[i] != 0 and flag == False):
            sign = -sign
            flag = true
            dprint("i: ", i, "change sign")

    dprint(" degree : ", abs(int(count/2)))
    return abs(int(count/2))

################################################################################
#RECURRENCE RELATIONS
################################################################################
def beta(num, symbol):
    if num < 0:
        return expand(-beta(-num,symbol))
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
        ret = (symbol*tau(num - 1, symbol) - tau(num - 2, symbol))
        return ret

def large_b(m,symbol1,symbol2):
    if np.any(m == 0):
        return 0
    dprint("large_b, input :", m)
    expr = 1
    multiplicand = 1
    for i in range(np.size(m)):
        if i%2 == 0:
            multiplicand = beta(m[i], symbol1)
            dprint("i:",i,", beta :", multiplicand)
            expr = expr * multiplicand
        else:
            multiplicand = beta(m[i], symbol2)
            dprint("i:",i,", beta :", multiplicand)
            expr = expr * multiplicand
    return expand(expr)

################################################################################
#MAIN LOGIC
################################################################################

def main():
    # input vector m and set dim
    print("\n type the value of m ")
    line = input(" ").split(' ')
    print("\n m : ", line)

    if len(line)%2 != 0:
        dim = len(line) + 1
    else:
        dim = len(line)

    r = (int)(dim/2)

    filename = "trace"

    m = np.zeros(dim, dtype=int)
    for i in range(len(line)):
        m[i] = (int)(line[i])
        filename = filename + "_" + str(line[i])

    filename = filename + ".txt"

    stdout = sys.stdout
    f = open(filename,'w', -1, "utf-8")
    sys.stdout = f
    
    print("m : ", line)

     # initiation
    init_printing(order='rev-lex')
    x,y,z = symbols('x,y,z')

    A=Matrix([[2,1],[0,1/2]])                                                       
    B=Matrix([[3,0],[1,1/3]]) 
    C=Matrix([[1,0],[0,1]])
    for i in range(dim):
        if i%2 == 0:
            C = C*A**(m[i])
        else:
            C = C*B**(m[i])

    # Compute expr1
    expr1 = expand((1/2) * tau(r,z) * large_b(m, x, y))
    
    # Compute expr2
    gen = eps_maker(dim)
    alt_count = 0
    expr2 = 0
    summand = 0
    deg = 0
    largeb = 0
    for mu in gen:
        if is_alternating(mu) == True:
            alt_count = alt_count+1
            print(alt_count,". mu:",  mu)
            deg = degree(mu)
            print(" degree of mu : ", deg,"\n")
            largeb = large_b(m + mu, x, y)
            if largeb == 0:
                summand = 0
                print( "-----------------------------------------------------------------------")
                continue
            else:
                summand = expand(((-1)**(r - deg)) * tau(deg,z) * largeb)
            dprint( " large_b : ", largeb )
            dprint( " sgn : ", ((-1)**(r - deg)) )
            dprint( " tau(deg) : ", tau(deg,z) )
            dprint( " summand : ", summand)
            expr2 = expr2 + summand       
            print( "-----------------------------------------------------------------------")
    
    if alt_count != (2**dim)-1:
        print("ERROR!!!")
    print(" expr1 : ", expr1)
    print(" expr2 : ", expr2)
    print( "-----------------------------------------------------------------------")
    expr = expand(expr1 + (1/2)*expr2)
    print(" RESULT is : ") 
    print(expr)
    print( "-----------------------------------------------------------------------")
    print( " Calculation Example ")
    print( " A : ", A)
    print( " B : ", B)
    print(" trace(C) : ", trace(C))
    print(" expr.sub : ", expr.subs({x:trace(A),y:trace(B),z:trace(A*B)}))
    
    f.close()
    sys.stdout = stdout
    print("")
    print(" RESULT is : \n ") 
    pretty_print(expr, order='rev-lex')
    print("\n Read " + filename + " for detailed calculation. ")


def tr(m):
    # input vector m and set dim
    x, y, z = symbols('x,y,z')
    dim = np.size(m)

    r = (int)((dim + 1)/ 2)

    print("size :" , dim, "r : ", r)
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
            deg = degree(mu)
            largeb = large_b(m + mu, x, y)
            if largeb == 0:
                summand = 0
                continue
            else:
                summand = expand(((-1) ** (r - deg)) * tau(deg, z) * largeb)
            dprint(" large_b : ", largeb)
            dprint(" sgn : ", ((-1) ** (r - deg)))
            dprint(" tau(deg) : ", tau(deg, z))
            dprint(" summand : ", summand)
            expr2 = expr2 + summand

    if alt_count != (2 ** dim) - 1:
        print("ERROR!!!")
    expr = expand(expr1 + (1 / 2) * expr2)

    return expr
################################################################################

if __name__ == "__main__":
    # initiation
    init_printing(order='rev-lex')
    A = Matrix([[2, 1], [0, 1 / 2]])
    B = Matrix([[3, 0], [1, 1 / 3]])
    C = Matrix([[1, 0], [0, 1]])

    expr1 = tr([3,9,6,18])
    expr2 = tr([3,9])

    print("expr1(f) :\n", expr1)
    print("\nexpr2(g) :\n", expr2)

    q,r = div(expr1,expr2)

    if (r == 0):
        print("\ndivided!!!, \n q :", q)
    else:
        print("\nNot divided!!")
