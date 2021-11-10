import numpy as np

#-----------------------------------------------------------------
# classes and functions for math
#-----------------------------------------------------------------

def poly(x,coeffs):
    '''return values for polynomial
    x is either a value or numpy array, coeffs=[c0,c1,c2...cn],
    return c0 + c1*x + c2*x^2 + ... + cn*x^n,
    '''
    # c0 + c1*x + c2*x^2 + ... + cn*x^n =
    # ((C[n]*x + C[n-1])*x + C[n-2])*x ... +C[0]
    # much faster than normal method
    l = len(coeffs)
    res = coeffs[l-1]
    for i in np.arange(l-1):
        res = res * x + coeffs[l-i-2]
    return res
    # res = 0.0
    # for i in np.arange(len(coeffs)):
    #     res += coeffs[i]*np.power(x,i)
    # return res

def chebyshev_poly(x,coeffs):
    """ return values for Fist Kind Chebyshev Polynomial
    x is either a value or numpy array, coeffs=[c0,c1,c2...cn],
    return c0*P[0] + c1*P[1] + c2*P[2] + ... + cn*P[n],
    where ord[i] is defined as:
        P[0] = 1
        P[1] = x
        P[2] = 2*x^2 - 1
        ......
        P[i] = 2*x*P[i-1] - P[i-2]
    """
    # append 0 for order 1 if only order 0
    if len(coeffs)==1:
        coeffs = np.array(coeffs)
        coeffs = np.append(coeffs,0)

    P   = [1.0,x]
    res = 0.0

    while(len(P)<len(coeffs)):
        P.append(2*x*P[-1]-P[-2])

    for i in np.arange(len(coeffs)):
        res += coeffs[i]*P[i]
    return res

def chebyshev2_poly(x,coeffs):
    """ return values for Second Kind Chebyshev Polynomial
    x is either a value or numpy array, coeffs=[c0,c1,c2...cn],
    return c0*P[0] + c1*P[1] + c2*P[2] + ... + cn*P[n],
    where ord[i] is defined as:
        P[0] = 1
        P[1] = 2x
        P[2] = 2*x^2 - 1
        ......
        P[i] = 2*x*P[i-1] - P[i-2]
    """
    # append 0 for order 1 if only order 0
    if len(coeffs)==1:
        coeffs = np.array(coeffs)
        coeffs = np.append(coeffs,0)

    P   = [1.0,2.0*x]
    res = 0.0

    while(len(P)<len(coeffs)):
        P.append(2*x*P[-1]-P[-2])

    for i in np.arange(len(coeffs)):
        res += coeffs[i]*P[i]
    return res


def legendre_poly(x,coeffs):
    """ return values for Legendre Polynomial
    x is either a value or numpy array, coeffs=[c0,c1,c2...cn],
    return c0*P[0] + c1*P[1] + c2*P[2] + ... + cn*P[n],
    where P[i] is defined as:
        P[0] = 1
        P[1] = x
        P[2] = 1/2*(3*x^2 - 1)
        ......
        P[i] = ( (2*i-1)*x*P[i-1] - (i-1)*P[i-2] )/i
    """
    # append 0 for order 1 if only order 0
    if len(coeffs)==1:
        coeffs = np.array(coeffs)
        coeffs = np.append(coeffs,0)

    P   = [1.0,x]
    res = 0.0

    i = 2
    while(len(P)<len(coeffs)):
        P.append(((2*i-1)*x*P[-1]-(i-1)*P[-2])/i)
        i += 1

    for i in np.arange(len(coeffs)):
        res += coeffs[i]*P[i]
    return res


