nele = 256
field_generator = 0x11d
gf_log = [0] * nele
gf_exp = [0] * nele * 2
code_generator = []
nsym = 4

def gf_init_table():
    global filed_generator
    global gf_log
    global gf_exp
    for i in range(nele - 1):
        if i == 0:
            tmp = 1
        else:
            tmp = tmp << 1
            if tmp >= nele:
                tmp = (tmp - nele) ^ (field_generator & (nele-1))
                
        gf_exp[i] = tmp
        gf_exp[i+nele-1] = tmp
        gf_log[tmp] = i
        
def gf_div(x,y):
    if y == 0:
        raise ZeroDivisionError()
    if x == 0:
        return 0
    return gf_exp[((gf_log[x] + nele - 1) - gf_log[y]) % (nele - 1)]

def gf_mul(x,y):
    if x == 0 or y == 0:
        return 0
    return gf_exp[gf_log[x] + gf_log[y]]

def gf_pow(x, power):
    if x == 0:
        return 0
    return gf_exp[(gf_log[x] * power) % (nele-1)]

def gf_inv(x):
    return gf_div(1, x)

def gf_poly_scale(p, x):
    if x == 1:
        return p
    elif x == 0:
        return [0]
    
    r = [0] * len(p)
    for i in range(len(p)):
        r[i] = gf_mul(p[i], x)
        
    return r

def gf_poly_add(p, q):
    r = [0] * max(len(p), len(q))
    
    for i in range(len(p)):
        r[i+len(r)-len(p)] = p[i]
        
    for i in range(len(q)):
        r[i+len(r)-len(q)] ^= q[i]
    return r
        
def gf_poly_mul(p, q):
    r = [0] * (len(p) + len(q) - 1)
    for j in range(len(q)):
        for i in range(len(p)):
            r[i+j] ^= gf_mul(p[i], q[j])
    return r

def gf_poly_eval(poly, x):
    y = poly[0]
    for i in range(1, len(poly)):
        y = gf_mul(y, x) ^ poly[i]
        
    return y

def gf_poly_div(dividend, divisor):
    normalizer = divisor[0]
    tmp_dividend = list(dividend) # copy of dividend
    tmp_divisor = gf_poly_scale(divisor, gf_inv(normalizer))
    for i in range(len(dividend) - len(divisor) + 1):
        if tmp_dividend[i] == 0:
            continue
        for j in range(1, len(tmp_divisor)):
            if tmp_divisor[j] == 0:
                continue
            
            tmp_dividend[i+j] ^= gf_mul(tmp_dividend[i], tmp_divisor[j])
            
    separator = -(len(tmp_divisor)-1)
    return gf_poly_scale(tmp_dividend[:separator], gf_inv(normalizer)), tmp_dividend[separator:]


def rs_init_code_generator():
    global code_generator
    global nsym
    code_generator = rs_generator_poly(nsym)

def rs_generator_poly(nsym):
    g = [1]
    for i in range(nsym):
        g = gf_poly_mul(g, [1, gf_pow(2, i)])
    return g


def rs_encode_msg(msg_in):
    global code_generator
    _, remainder = gf_poly_div(msg_in + [0] * (len(code_generator)-1), code_generator)
    msg_out = msg_in + remainder
    return msg_out


def rs_calc_syndromes(msg):
    synd = [0] * (len(code_generator) - 1)
    for i in range(len(synd)):
        synd[len(synd) - i - 1] = gf_poly_eval(msg, gf_pow(2,i))

    return synd

def rs_check_syndromes(synd):
    return max(synd) == 0

def rs_find_error_locator_and_evaluator(synd):
    v_pre = [0] * (nsym + 1)
    v_pre[0] = 1
    v_cur = list(synd)
    x_pre = [0]
    x_cur = [1]
    while True:
        if len(v_cur) <= (nsym/2):
            break
        
        q,r = gf_poly_div(v_pre, v_cur)
        v_pre = v_cur
        v_cur = r
        x_tmp = x_pre
        x_pre = x_cur
        x_cur = gf_poly_add(x_tmp, gf_poly_mul(q, x_cur))
    return x_cur, v_cur

def rs_correct_errata(r, locator, evaluator):
    err_loc = rs_find_errors(locator, len(r))
    msg_in = list(r)
    dlocator = [locator[i] for i in range(1,len(locator),2)]
    for e in err_loc:
        x = gf_exp[e]
        x_inv = gf_inv(x)
        mag = gf_mul(x, gf_div(gf_poly_eval(evaluator, x_inv), gf_poly_eval(dlocator, x_inv)))
        msg_in[len(r)-e-1] = mag ^ r[len(r)-e-1]

    return msg_in

def rs_find_errors(locator, nmsg):
    ret = []
    a = gf_inv(gf_exp[nele - nmsg])
    tmp_locator = [0] * len(locator)
    for i in range(len(locator)):
        tmp_locator[i] = gf_mul(locator[i], gf_pow(a, i))
    
    for i in range(nmsg):
        res = reduce(lambda x,y: x^y, tmp_locator)
        for j in range(1, len(tmp_locator)):
            index = len(tmp_locator) - j - 1
            tmp_locator[index] = gf_mul(tmp_locator[index], gf_exp[j])
        if res == 0:
            ret.append(nmsg - i - 1)

    return ret

if __name__ == "__main__":
    print("")
    gf_init_table()
    rs_init_code_generator()
    encoded_data = rs_encode_msg([1,2,3,4,5,6,7,8,9,10,11])
    correct_synd = rs_calc_syndromes(encoded_data)
    if not rs_check_syndromes(correct_synd):
        raise RuntimeError("syndrome should be an all zero array when the data is correct")
    
    err = [0,0,0,0,0,13,0,0,0,0,0,0,2,0,0]
    r = gf_poly_add(encoded_data, err)
    synd = rs_calc_syndromes(r)
    locator, evaluator = rs_find_error_locator_and_evaluator(synd)
    
    msg = rs_correct_errata(r, locator, evaluator)
    print("")
    print("result",msg)
    print("expect",encoded_data)
    # test for poly div

    
