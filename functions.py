from math import sqrt,log

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

def calculate(x,find = 0.0):
	if x < 0:
		z = (-x*H(x)**5)/G(x)**5
		rho = -x*G(x)*H(x**6)*P(x**3)/P(x)
		prod = 1.0
		for i in range(1,times,1):
			prod = prod*((1.0-x**(6.0*i-4.0))*(1.0-x**(6.0*i-3.0))**2*(1.0-x**(6.0*i-2.0)))/((1.0-x**(6.0*i-5.0))*(1.0-x**(6.0*i-1.0))*(1.0-x**(6.0*i))**2)
		btp = H(x)**3*Q(x**5)**2/G(x)**2*prod
	else:
		z = G(x)**5/((x)*H(x)**5)
		rho23 = x**2*H(x)*Q(x)*H(x**9)*Q(x**9)/Q(x**3)**2
		rho1 = H(x)*Q(x)*(G(x)*Q(x)+x**2*H(x**9)*Q(x**9))/Q(x**3)**2
		rho = (rho1 + rho23*2.0)/3.0
		prod = 1.0
		for i in range(1,times,1):
			prod = prod*((1.0-x**(3.0*i-2.0))*(1.0-x**(3.0*i-1.0)))/(1.0-x**(3.0*i))**2
		btp = x**(-1.0/3.0)*G(x)**3*Q(x**5)**2/H(x)**2*prod
	return log(z)-find, -log(btp),rho
