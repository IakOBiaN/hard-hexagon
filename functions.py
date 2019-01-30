from math import sqrt

#iterations count
times = 1000

def G(x,times=times):
    prod = 1.
    for i in range(1,times):
        prod = prod*1./((1.-x**(5.*i-4.))*(1.-x**(5.*i-1.)))
    return prod

def H(x,times=times):
    prod = 1.
    for i in range(1,times):
        prod = prod*1./((1.-x**(5.*i-3.))*(1.-x**(5.*i-2.)))
    return prod

def P(x,times=times):
    prod = 1.
    for i in range(1,times):
        prod = prod*(1.-x**(2.*i-1.))
    return prod

def Q(x,times=times):
    prod = 1.
    for i in range(1,times,1):
        prod = prod*(1.-x**i)
    return prod
