from sympy import * 
import numpy as np

################################################################################
# FILE LOGGING
################################################################################
''' print into the logfile'''
def fprint( *args):
    print(*args, file=fileinfo)
    return

################################################################################
# EPSILON RELATED
################################################################################

def eps_maker(dim):
    ''' Using generator, it makes all possible epsilons '''
    fprint("start eps_generator, dim: ",dim)
    print( "-----------------------------------------------------------------------")
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
    fprint("finish eps_generator")
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
    eps = cut_out(eps)
    
    # Roll left until first entry becomes 3
    eps = roll_until(eps, 3) 

    for i in range(size):
        # start of subinterval
        fprint("eps[i] : ", eps[i],", ", flag, sgn_flag)
        if eps[i] == 0 and flag == False:
            sgn = 1
            subcount = 1
            flag = True
            sgn_flag = False
                
        # end of subinterval, or inside [3,3,...,3]
        elif eps[i] == 3 :
            count = count + abs(subcount)
            fprint("count : ", count)
            subcount = 0
            flag = False
        # zero, inside a subinterval
        # since sgn can change after 0, reset sgn_flag
        elif eps[i] == 0 and flag == True:
            if sgn_flag == True:
                sgn = -sgn
                fprint("change sgn")
                sgn_flag = False
            subcount = subcount + sgn
        # nonzero, inside a subinterval
        elif eps[i] != 0 and flag == True:
            sgn_flag = True

        # last entry
        if i == size - 1:
                count = count + abs(subcount)
                fprint("count : ", count)

    return (int)(1/2 * count)

################################################################################
# OPERATIONS ON EPSILON
################################################################################


def roll_until(eps, num):
    ''' roll(-1) [eps] until first entry becomes [num] '''
    size_eps = np.size(eps) 

    if np.all(eps != num) == True:
        fprint("no entry is ", num," stop rolling. ")
        return eps

    while (eps[0] != num):
        eps = np.roll(eps,-1) 
    fprint('rolling: ',eps)
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
                            eps[i-j-1] = 3
                            fprint("1 middle cut out : ", eps, "even leng: ", even_leng, "i: ", i)
                            count = 0
                    else:
                        count = 0
                elif i == size - 1:
                    if count == even_leng:
                        for j in range(count):
                            eps[i-j] = 3
                            fprint("2 middle cut out : ", eps, "even leng: ", even_leng, "i: ", i)
                            count = 0
        even_leng = even_leng - 2
    fprint("cut out : ", eps)
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
            fprint("i: ",i,"count: ",count)
        elif (sub_eps[i] != 0 and flag == False):
            sign = -sign
            flag = true
            fprint("i: ", i, "change sign")

    fprint(" degree : ", abs(int(count/2)))
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
    fprint("large_b, input :", m)
    expr = 1
    multiplicand = 1
    for i in range(np.size(m)):
        if i%2 == 0:
            multiplicand = beta(m[i], symbol1)
            fprint("i:",i,", beta :", multiplicand)
            expr = expr * multiplicand
        else:
            multiplicand = beta(m[i], symbol2)
            fprint("i:",i,", beta :", multiplicand)
            expr = expr * multiplicand
    return expand(expr)

################################################################################
#MAIN LOGIC
################################################################################

def main():
    # input vector m and set dim
    print("\n type the value of m ")
    line = input(" ").strip().split(' ')
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

    global fileinfo
    fileinfo = open(filename,'w', -1, "utf-8")

    fprint("m : ", line)

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
            fprint(alt_count,". mu:",  mu)
            deg = degree(mu)
            fprint(" degree of mu : ", deg,"\n")
            largeb = large_b(m + mu, x, y)
            if largeb == 0:
                summand = 0
                fprint( "-----------------------------------------------------------------------")
                continue
            else:
                summand = expand(((-1)**(r - deg)) * tau(deg,z) * largeb)
            fprint(" large_b : ", largeb )
            fprint(" sgn : ", ((-1)**(r - deg)) )
            fprint(" tau(deg) : ", tau(deg,z) )
            fprint(" summand : ", summand)
            expr2 = expr2 + summand       
            fprint( "-----------------------------------------------------------------------")
    
    if alt_count != (2**dim)-1:
        print("ERROR!!!")
    fprint(" expr1 : ", expr1)
    fprint(" expr2 : ", expr2)
    fprint( "-----------------------------------------------------------------------")
    expr = expand(expr1 + (1/2)*expr2)
    fprint(" RESULT is : ")
    fprint(expr)
    fprint( "-----------------------------------------------------------------------")
    fprint( " Calculation Example ")
    fprint( " A : ", A)
    fprint( " B : ", B)
    fprint(" trace(C) : ", trace(C))
    fprint(" expr.sub : ", expr.subs({x:trace(A),y:trace(B),z:trace(A*B)}))
    
    print("")
    print(" RESULT is : \n ") 
    pretty_print(expr, order='rev-lex')
    print("\n Read " + filename + " for detailed calculation. ")

################################################################################

if __name__ == "__main__":
    main()
