# -*- coding: utf-8 -*-
"""
Created on Tue Jan 29 23:14:17 2013

@authors: David Kingman, Valentino Constantinou, Taylor Stevens
"""

#NOTE: There are currently some errors in the code (linear algebra). 

import math
import copy
from numarray import argmax
from math import cos
from math import atan
from math import sin
from math import pi
from math import 
from numpy import array

def eulersConstant():
    x = s = 1.0
    for i in range(5000):
        if ((i % 100) == 0):
            print "x = " + str(x) + " s = " + str(s)
        x = x + 1.0
        s = s + 1.0 / x
    print "x = " + str(x) + " s = " + str(s)


eulersConstant()
print 'If the program ran for a week it might reach the machine maximum. The series has logarithmic growth so it will grow very slowly initialy\n'


def absError(real, approx):
    return approx - real


def relError(real, approx):
    return approx / real


def my_exp(x, precision=1e-6, max_steps=25):
    if x == 0:
        return 1.0
    elif x > 0:
        return 1.0 / my_exp(-x, precision, max_steps)
    else:
        t = s = 1.0
        for k in range(1, max_steps):
            t = t * x / k
            s = s + t
            print s
            if abs(t) < precision:
                print 'Approx = ', s, ' Exact = ', math.exp(x)
                print 'n= ', k, ' Abs Error= ', absError(math.exp(x), s), ' Rel Error= ', relError(math.exp(x), s)
                return s
        raise ArithmeticError, 'no convergence'


my_exp(0)
my_exp(1)
my_exp(-1)
my_exp(.5)
my_exp(-.123)


def func_y1(x):
    return x * (x * (x - 1) - 2) + 1


def bisection(f, a, b, ap=1e-6, max_steps=50):
    n = 1.0
    while n <= max_steps:
        c = (a + b) / 2.0
        #print c
        if f(c) == 0 or (b - a) / 2 < ap:
            return c
        n += 1
        if (f(c) > 0 and f(a) > 0) or (f(c) < 0 and f(a) < 0):
            a = c
        else:
            b = c
    if n == max_steps:
        raise ArithmeticError, 'Method failed'


print 'Result for x^3-x^2-2x+1 in interval (0,1): ', bisection(func_y1, 0.0, 1.0)
print 'Result for x^3-x^2-2x+1 in interval (1,2): ', bisection(func_y1, 1.0, 2.0)
print


def func_y2(x):
    return 9.0 * x ** 4 + 18 * x ** 3 + 38 * x ** 2 - 54 * x + 14


print 'Result for 9x^4+18x^3+38x^2-54x+14 in interval (0,1): ', bisection(func_y2, 0.0, 1.0)


def bisection_2(ap=1e-6, max_steps=50):
    f = lambda x: x * (x * x + 3) - 1
    g = lambda x: x ** 3 * math.sin(x)
    h = lambda x: x + 10 - math.cosh(50 / x)
    n = 1.0
    a = 0.0
    b = 1.0
    while n <= max_steps:
        c1 = (a + b) / 2.0
        if f(c1) == 0 or (b - a) / 2 < ap:
            print 'Root for f(x): ', c1
            break
        n += 1
        if (f(c1) > 0 and f(a) > 0) or (f(c1) < 0 and f(a) < 0):
            a = c1
        else:
            b = c1
    if n == max_steps:
        raise ArithmeticError, 'Method failed'
    a = .5
    b = 2.0
    while n <= max_steps:
        c2 = (a + b) / 2.0
        if f(c2) == 0 or (b - a) / 2 < ap:
            print 'Root for g(x): ', c2
            break
        n += 1
        if (g(c2) > 0 and g(a) > 0) or (g(c2) < 0 and g(a) < 0):
            a = c2
        else:
            b = c2
    if n == max_steps:
        raise ArithmeticError, 'Method failed'
    a = 120.0
    b = 130.0
    while n <= max_steps:
        c3 = (a + b) / 2.0
        if h(c3) == 0 or (b - a) / 2 < ap:
            print 'Root for h(x): ', c3
            break
        n += 1
        if (h(c3) > 0 and h(a) > 0) or (h(c3) < 0 and h(a) < 0):
            a = c3
        else:
            b = c3
    if n == max_steps:
        raise ArithmeticError, 'Method failed'


bisection_2()

def newton(f,f_prime,x,tolerance,precision,display,steps):
    fx=f(x)
    for i in range(steps):
        fp=f_prime(x)
        if abs(fp) < precision:
            display=0
            print 'small derivative'
            #return [steps,x,fx]
            break
        d=fx/fp
        x=x-d
        fx = f(x)
        if display:
                print "n = %i, x = %f, and xn=%f"% (i + 1, x, fx)
                #print [steps,x,fx]
        if abs(d) < tolerance:
            print 'convergence'
            break
    return [steps, x, fx]


def function(x):
    return (tan(x))-x

def functionprime(x):
    return (tan(x)**2)

#def functionprimeprime(x):
    #return 2*(tan(x))*((sec(x)**2))

print newton(function,functionprime,7,1*e**-6,1*e**-6,1,25)

#None. No root. This answer makes sense when you look at the plot.

def function2(x):
    return (e**(x))-((x+9)**(0.5))

def function2prime(x):
    return (e**(x))-(1/(2*(x+9)**(0.5))) #may need to rewrite this function to minimize error.

#def function2primeprime(x):
    #return (e**(x)) - 1/(4*((x+9)**(3/2))) #may need to rewrite this function to minimize error.

print newton(function2,function2prime,2,1*e**-6,1*e**-6,1,25)

def newtonmod(f,f_prime,x,tolerance,precision,display,steps):
    fx=f(x)
    for i in range(steps):
        fp=f_prime(x)
        if abs(fp) < precision:
            display=0
            print 'small derivative'
            #return [steps,x,fx]
            break
        d=fx/fp
        x=x-d
        if abs(f(x-d)) >= abs(f(x)):
            d=0.5*d
        else:
            d=d
        fx = f(x)
        if display:
                print "n = %i, x = %f, and xn=%f"% (i + 1, x, fx)
        #print [steps,x,fx]
        if abs(d) < tolerance:
            print 'convergence'
            break
    return [steps, x, fx]



def testfn(x):
    return sin(x)

def testfnprime(x):
    return cos(x)


print newton(testfn,testfnprime,1.2,1*e**-6,1*e**-3,1,25)
print newtonmod(testfn,testfnprime,1.2,1*e**-6,1*e**-3,1,25)

#yay it works.

def newtonaccel(f,f_prime,x,tolerance,precision,display,steps):
    fx=f(x)
    for i in range(steps):
        fp=f_prime(x)
        if abs(fp) < precision:
            display=0
            print 'small derivative'
            #return [steps,x,fx]
            break
        d=fx/fp
        x=x-(2*d)
        #if abs(f(x-d)) >= abs(f(x)):
         #   d=0.5*d
        #else:
         #   d=d
        fx = f(x)
        if display:
                print "n = %i, x = %f, and xn=%f"% (i + 1, x, fx)
        #print [steps,x,fx]
        if abs(d) < tolerance:
            print 'convergence'
            break
    return [steps, x, fx]


def doublerootfn(x):
    return (e**x)*((x-1)**2)

def doublerootfnprime(x):
    return (e**x)*((x**2)-1)

def secant(f, a, b, ap=1e-6, max_steps=25):
    fa = f(a)
    fb = f(b)
    if abs(fa) > abs(fb):
        temp = a
        a = b
        b = temp
        temp = fa
        fa = fb
        fb = temp
    for n in range(2, max_steps):
        print 'n=', n, ' a=', a, ' fa=', fa
        if abs(fa) > abs(fb):
            temp = a
            a = b
            b = temp
            temp = fa
            fa = fb
            fb = temp
        d = (b - a) / (fb - fa)
        b = a
        fb = fa
        d = d * fa
        if abs(d) < ap:
            print 'Convergence'
            return [n, a, fa]
        a = a - d
        fa = f(a)
        return [n, a, fa]


#secant method
def secant(f, a, b, precision, steps):
    fa = f(a)
    fb = f(b)
    if abs(fa) > abs(fb):
        temp = a
        a = b
        b = temp
        temp = fa
        fa = fb
        fb = temp
    print [0,a,fa]
    print [1,b,fb]
    for n in range(2, steps):
        if abs(fa) > abs(fb):
            temp = a
            a = b
            b = temp
            temp = fa
            fa = fb
            fb = temp
        d = (b - a) / (fb - fa)
        b = a
        fb = fa
        d = d * fa
        if abs(d) < precision:
            print 'Convergence'
            #print "" #so that when testing multiple functions the output is easier to read.
            break
        #else: #Do not need this.
            #print 'Does Not Converge'
            #print "" #so that when testing multiple functions the output is easier to read.
            #break
        a = a - d
        fa = f(a)
        print 'n=', n, ' a=', a, ' fa=', fa
    return [n, a, fa]


def secantfn(x):
    return (e**x)-3*(x**2)

print secant(secantfn,-0.5,2,1*e**-10,25)
#found root to be -0.45896257524. There are other roots at 0.91 and 3.733

def wilkinson(x):
    #n=20
    #for i in range(n):
        #return(sum(x-n))
        return (x-1)*(x-2)*(x-3)*(x-4)*(x-5)*(x-6)*(x-7)*(x-8)*(x-9)*(x-10)*(x-11)*(x-12)*(x-13)*(x-14)*(x-15)*(x-16)*(x-17)*(x-18)*(x-19)*(x-20)

print secant(wilkinson,20,21,1*e**-10,100)
#here, the sequence converges to 20

def wilkinsonmod(x):
    return ((x-1)*(x-2)*(x-3)*(x-4)*(x-5)*(x-6)*(x-7)*(x-8)*(x-9)*(x-10)*(x-11)*(x-12)*(x-13)*(x-14)*(x-15)*(x-16)*(x-17)*(x-18)*(x-19)*(x-20)) - ((10**-8)*(x**19))

print secant(wilkinsonmod,20,21,1*e**-10,100)
#here, the sequence does 20.240274, between [20,21]. 

f1 = lambda x: (1. / math.e ** (x * x)) - math.cos(x) - 1.
f1_p = lambda x: math.sin(x) - 2. * x * (1. / math.e ** (x * x))

print 'x=0 ', newton(f1, f1_p, 0.)
print 'x=1 ', newton(f1, f1_p, 1.)
print 'We see when we pick 3 and 4 as our x_0 we get the two correct roots in this range'
print 'x=4 ', newton(f1, f1_p, 4.)
print 'x=3 ', newton(f1, f1_p, 3.)


def newton2(f, f_prime, x, delta=1e-6, ap=1e-6, max_steps=25):
    fx = f(x)
    for n in range(1, max_steps):
        fp = f_prime(x)
        if abs(fp) < delta:
            print 'Small derivative'
            return [n, x, fx]
        d = 2 * fx / fp
        x = x - d
        fx = f(x)
        if abs(d) < ap:
            print 'Convergence'
            return [n, x, fx]


f2 = lambda x: (x - 1.) ** 2.
f2_p = lambda x: 2. * x - 2.
print newton2(f2, f2_p, 2.)
print newton(f2, f2_p, 2.)
print 'We can see it took twenty steps in the original Newtons method, but only 2 in the modified version'

f_oliver = lambda x: x * x - 1
f_oliver_p = lambda x: 2. * x
f_oliver_pp = lambda x: 2.


def oliver(f, f_prime, f_prime_prime, x, ap=1e-6, max_steps=25):
    fx = f(x)
    for n in range(1, max_steps):
        fp = f_prime(x)
        fpp = f_prime_prime(x)
        if abs(fp) < ap:
            print 'Small derivative'
            return [n, x, fx]
        fd = fx / fp - .5 * fpp / fp * (fx / fp) ** 2
        d = x - (x - fd)
        x = x - fd
        fx = f(x)
        if abs(d) < ap:
            print 'Convergence'
            return [n, x, fx]


print 'Our function f(x) = x^2-1 has a root at 1'
print 'We see that when we start at x_0 = 2 we get to 1 with a precision'
print 'of 1e-6 in six steps'
print oliver(f_oliver, f_oliver_p, f_oliver_pp, 2.)


def coef(n, xi, yi, ai):
    for i in range(0, n):
        ai[i] = yi[i]
    for j in range(1, n):
        for i in range(n - 1, j - 1, -1):
            ai[i] = (ai[i] - ai[i - 1]) / (xi[i] - xi[i - j])
    return ai


xi = [1., 2., 3., -4., 5]
yi = [2., 48., 272., 1182., 2262.]
n = len(yi)
ai = [0 for i in range(n)]

a = coef(n, xi, yi, ai)
print 'Coeffecients are: ', a


def eval(n, x, a, t):
    temp = a[n - 1]
    for i in range(n - 2, -1, -1):
        temp = temp * (t - x[i]) + a[i]
    return temp


t = -1.
result = eval(n, xi, a, t)
print 'Polynomial evaluated at ', t, ' is ', result

f3 = lambda x: math.exp(1) ** x
t = 1.
xi = [0, .2, .4, .6, .8, 1., 1.2, 1.4, 1.6, 1.8, 2.]
yi = [f3(i) for i in xi]
n = len(yi)
ai = [0 for i in range(n)]
a = coef(n, xi, yi, ai)
result = eval(n, xi, a, t)
for i in range(0, 100):
    t = i
    result = eval(n, xi, a, t)
    print 'e(', i, ') evaluated using polynomial: ', t, ' is ', result, ' Real: ', math.exp(i)


def coef(n, xi, yi, ai):
    for i in range(0, n):
        ai[i] = yi[i]
    for j in range(1, n):
        for i in range(n - 1, j - 1, -1):
            ai[i] = (ai[i] - ai[i - 1]) / (xi[i] - xi[i - j])
    return ai


def eval(n, x, a, t):
    temp = a[n - 1]
    for i in range(n - 2, -1, -1):
        temp = temp * (t - x[i]) + a[i]
    return temp


def func_1(x): return 1. / (x * x + 1.)


xi = np.arange(-5, 5.5, .5)
yi = [func_1(x) for x in xi]
n = len(yi)
ai = [0 for i in range(n)]
a = coef(n, xi, yi, ai)

xi_2 = np.arange(-5, 5.25, .25)
print "{:11s}".format(' '), "{:^12s}".format("p(x)"), "{:^14s}".format("f(x)")
for x in xi_2:
    approx = eval(n, xi, a, x)
    real = func_1(x)
    print 'x =', "{: 2.2f}".format(x), "{: 1.10f}".format(approx), "{: 1.10f}".format(real)
print "Notice the wild flucuations on the nodes that weren't part of the original 21 nodes"


def cheb_nodes(xi, a, b, n):
    return [.5 * (a + b) + .5 * (b - a) * math.cos((2. * i + 1.) / (2. * n + 2.) * math.pi) for i in range(0, n)]


chebnodes = cheb_nodes(xi, -5., 5., 21)
chebnodes.reverse()

yi = [func_1(x) for x in chebnodes]
n = len(yi)
ai = [0 for i in range(n)]
a = coef(n, chebnodes, yi, ai)
print "{:11s}".format(' '), "{:^12s}".format("p(x)"), "{:^14s}".format("f(x)")
for x in xi_2:
    approx = eval(n, chebnodes, a, x)
    real = func_1(x)
    print 'x =', "{: 2.2f}".format(x), "{: 1.10f}".format(approx), "{: 1.10f}".format(real)
print 'Notice that all the approximations are now much closer than equally spaced nodes'


def func_2(x): return abs(x)


xi = np.arange(-1, 1.2, .2)
yi = [func_2(x) for x in xi]
n = len(yi)
ai = [0 for i in range(n)]
a = coef(n, xi, yi, ai)
xi_2 = np.arange(-1, 1.05, .05)

print'\nEqually Spaced Nodes'
print "{:^14s}".format('|x| - p(x)')
for x in xi_2:
    approx = eval(n, xi, a, x)
    real = func_2(x)
    print "{: 1.12f}".format(real - approx)

chebnodes = cheb_nodes(xi, -1., 1., 11)
chebnodes.reverse()
yi = [func_2(x) for x in xi]
n = len(yi)
ai = [0 for i in range(n)]
a = coef(n, chebnodes, yi, ai)

print'\nChebyshev Nodes'
print "{:^14s}".format('|x| - p(x)')
for x in xi_2:
    approx = eval(n, chebnodes, a, x)
    real = func_2(x)
    print "{: 1.12f}".format(real - approx)




def createA(n): #function to create matrix A
    A = np.empty((n,n), order='F') #declare an empty array of nxn dimension
    for i in range(n):
        for j in range(n):
             if j>=1:
                 A[1,j]=A[j,1]=n**(-1)
             else:
                 A[i,j]=A[i-1,j]+A[i,j-1]
    return A

A = createA(25)
#print A

def det(a):
    n = len(a)
    for k in range(0,n-1):
        for i in range(k+1,n):
            if a[i,k] != 0.0:
                lam = a [i,k]/a[k,k]
                a[i,k+1:n] = a[i,k+1:n] - lam*a[k,k+1:n]
                a[i,k] = lam
    det = product(diagonal(a)) #a triangular matrix det is the product of the diagonal
    return det #toggle these two return functions to give either the determinant or triangularly reduced A.
    #return a

detA = det(A)
print detA

#Here we can clearly see that the determinant is 0 or undefined. Looking at the reduced A,
#we see that there are zeros in the diagonal and entire blocks of zeros in the matrix... This won't be easily solveable.

B = createA(6)
print B

#same happens on smaller scale. 

def createA(n): #function to create matrix A
    A = np.empty((n,n), order='F') #declare an empty array of nxn dimension
    for i in range(n):
        for j in range(n):
             if j>=i:
                 A[i,j]=1
             else:
                 A[i,j]=-1
    return A

A = createA(25)
print A

def det(a):
    n = len(a)
    for k in range(0,n-1):
        for i in range(k+1,n):
            if a[i,k] != 0.0:
                lam = a [i,k]/a[k,k]
                a[i,k+1:n] = a[i,k+1:n] - lam*a[k,k+1:n]
                a[i,k] = lam
    det = product(diagonal(a)) #a triangular matrix det is the product of the diagonal
    return det #toggle these two return functions to give either the determinant or triangularly reduced A.
    #return a

detA = det(A)
print detA

#here we obtain a determinant of 1677216.

def createA(n): #function to create matrix A
    A = np.empty((n,n), order='F') #declare an empty array of nxn dimension
    for i in range(n):
        for j in range(n):
             A[i,j]= abs(i-j)
    return A

A = createA(25)
#print A

def det(a):
    n = len(a)
    for k in range(0,n-1):
        for i in range(k+1,n):
            if a[i,k] != 0.0:
                lam = a [i,k]/a[k,k]
                a[i,k+1:n] = a[i,k+1:n] - lam*a[k,k+1:n]
                a[i,k] = lam
    det = product(diagonal(a))
    #return A #toggle these two return functions to give either the determinant or triangularly reduced A.
    return det

detA = det(A)
print detA

#Here we can clearly see that the determinant is 0 or undefined. Looking at the reduced A,
#we see that there are zeros in the diagonal.

def createA(n): #function to create matrix A
    A = np.empty((n,n), order='F') #declare an empty array of nxn dimension
    for i in range(n):
        for j in range(n):
             A[i,j]= -1+2*argmax([i,j])
    return A

A = createA(30)
print A

def createB(n): #function to create matrix B
    B = np.empty((n,1), order='F') #declare an empty array of nxn dimension
    for i in range(n):
        for j in range(n):
             B[i] = 1
    return B

B = createB(30)
print B

def gaussA(n, a):
    s = [0 for i in range(n)]
    x = [0 for i in range(n)]
    l = [i for i in range(n)]
    for i in range(0, n):
        l[i] = i
        smax = 0
        for j in range(0, n):
            #print abs(a[i,j])
            smax = max(smax, abs(a[i, j]))
        s[i] = smax
        #print s
        #print
    for k in range(0, n - 1):
        rmax = 0
        for i in range(k, n):
            #print k/s[l[i]]
            r = abs(a[l[i], l[k]] / s[l[i]])
            #print r>rmax
            #print 'j=',j,' i=',i, 'r=',r,'rmax=',rmax
            if r > rmax:
                rmax = r
                j = i

        #print k
        temp = l[k]
        l[k] = l[j]
        l[j] = temp
        #print l
        #print '\n'
        for i in range(n - 1, k, -1):
            xmult = a[l[i], k] / a[l[k], k]
            #print xmult
            a[l[i], k] = xmult
            for j in range(k + 1, n):
                a[l[i], j] = a[l[i], j] - xmult * a[l[k], j]
    return a

def gaussl(n, a):
    s = [0 for i in range(n)]
    x = [0 for i in range(n)]
    l = [i for i in range(n)]
    for i in range(0, n):
        l[i] = i
        smax = 0
        for j in range(0, n):
            #print abs(a[i,j])
            smax = max(smax, abs(a[i, j]))
        s[i] = smax
        #print s
        #print
    for k in range(0, n - 1):
        rmax = 0
        for i in range(k, n):
            #print k/s[l[i]]
            r = abs(a[l[i], l[k]] / s[l[i]])
            #print r>rmax
            #print 'j=',j,' i=',i, 'r=',r,'rmax=',rmax
            if r > rmax:
                rmax = r
                j = i

        #print k
        temp = l[k]
        l[k] = l[j]
        l[j] = temp
        #print l
        #print '\n'
        for i in range(n - 1, k, -1):
            xmult = a[l[i], k] / a[l[k], k]
            #print xmult
            a[l[i], k] = xmult
            for j in range(k + 1, n):
                a[l[i], j] = a[l[i], j] - xmult * a[l[k], j]
    return l

def solve(n, a, l, b): #takes a matrix A, b, l (the scale factors), of size n and solves.
    x = [0 for i in range(n)]
    for k in range(0, n - 1):
        for i in range(k + 1, n):
            b[l[i]] = b[l[i]] - a[l[i], k] * b[l[k]]
            #print n,x
    x[n - 1] = b[l[n - 1]] / a[l[n - 1], n - 1]
    for i in range(n - 1, -1, -1):
        s = b[l[i]]
        for j in range(i + 1, n):
            s = s - a[l[i], j] * x[j]
        x[i] = s / a[l[i], i]
    return x

G = gaussA(30,A)
L = gaussl(30,A)
X = np.asmatrix(solve(30,A,L,B), dtype=np.int8)
print X

#and we obtain an end result.


def createHilbert(n): #define a function for creating Hilbert matrices.
    H = np.empty((n, n), order='F') #declare an empty array of nxn dimension.
    for j in range(n):
        for i in range(n):
            H[i,j] = 1. / (i + j + 1) #formula for constructing matrix values.
    return H #return the matrix A

H3 = createHilbert(3) #create 3x3 Hilbert matrix
print H3

H8 = createHilbert(8) #create 8x8 Hilbert matrix
#print H8

H13 = createHilbert(13) #create 13x13 Hilbert matrix
#print H13

def createB(n): #function to create matrix b
    H = np.empty((n, n), order='F') #declare an empty array of nxn dimension.
    for j in range(n):
        for i in range(n):
            H[i,j] = 1. / (i + j + 1) #formula for constructing matrix values.
    B = np.empty((n,1), order='F') #declare an empty array of nxn dimension
    for i in range(n):
        B[i]=0
        for j in range(n):
            B[i]=B[i]+H[i,j] #formula for constructing matrix values.
    return B

B3 = createB(3)
print B3 #create 3x1 matrix of B's using the formula.

B8 = createB(8)
#print B8 #create 8x1 matrix of B's using the formula.

B13 = createB(13)
#print B13 #create 13x1 matrix of B's using the formula.

def gaussA(n, a):
    s = [0 for i in range(n)]
    x = [0 for i in range(n)]
    l = [i for i in range(n)]
    for i in range(0, n):
        l[i] = i
        smax = 0
        for j in range(0, n):
            #print abs(a[i,j])
            smax = max(smax, abs(a[i, j]))
        s[i] = smax
        #print s
        #print
    for k in range(0, n - 1):
        rmax = 0
        for i in range(k, n):
            #print k/s[l[i]]
            r = abs(a[l[i], l[k]] / s[l[i]])
            #print r>rmax
            #print 'j=',j,' i=',i, 'r=',r,'rmax=',rmax
            if r > rmax:
                rmax = r
                j = i

        #print k
        temp = l[k]
        l[k] = l[j]
        l[j] = temp
        #print l
        #print '\n'
        for i in range(n - 1, k, -1):
            xmult = a[l[i], k] / a[l[k], k]
            #print xmult
            a[l[i], k] = xmult
            for j in range(k + 1, n):
                a[l[i], j] = a[l[i], j] - xmult * a[l[k], j]
    return a

def gaussl(n, a):
    s = [0 for i in range(n)]
    x = [0 for i in range(n)]
    l = [i for i in range(n)]
    for i in range(0, n):
        l[i] = i
        smax = 0
        for j in range(0, n):
            #print abs(a[i,j])
            smax = max(smax, abs(a[i, j]))
        s[i] = smax
        #print s
        #print
    for k in range(0, n - 1):
        rmax = 0
        for i in range(k, n):
            #print k/s[l[i]]
            r = abs(a[l[i], l[k]] / s[l[i]])
            #print r>rmax
            #print 'j=',j,' i=',i, 'r=',r,'rmax=',rmax
            if r > rmax:
                rmax = r
                j = i

        #print k
        temp = l[k]
        l[k] = l[j]
        l[j] = temp
        #print l
        #print '\n'
        for i in range(n - 1, k, -1):
            xmult = a[l[i], k] / a[l[k], k]
            #print xmult
            a[l[i], k] = xmult
            for j in range(k + 1, n):
                a[l[i], j] = a[l[i], j] - xmult * a[l[k], j]
    return l

G3 = gaussA(3, H3) #Gaussian elimination with scaled pivoting is carried out on our Hilbert matrices.
#print G3

L3 = gaussl(3,H3)
L3 = np.array(L3)

G8 = gaussA(8, H8) #Gaussian elimination with scaled pivoting is carried out on our Hilbert matrices.
#print G8
G13 = gaussA(13, H13) #Gaussian elimination with scaled pivoting is carried out on our Hilbert matrices.
#print G13

def solve(n, a, l, b): #takes a matrix A, b, l (the scale factors), of size n and solves.
    x = [0 for i in range(n)]
    for k in range(0, n - 1):
        for i in range(k + 1, n):
            b[l[i]] = b[l[i]] - a[l[i], k] * b[l[k]]
            #print n,x
    x[n-1] = b[l[n - 1]] / a[l[n - 1], n - 1]
    for i in range(n - 1, -1, -1):
        s = b[l[i]]
        for j in range(i + 1, n):
            s = s - a[l[i], j] * x[j]
            x[i] = s / a[l[i], i]
    return x

#Now we need to solve for X3, X8, and X13 respectively.

#X3 = solve(3,H3,L3,B3)
#print X3
X3 = np.asmatrix(solve(3,H3,L3,B3), dtype=np.uint8) #need to do this to display properly.
print X3

#ones all around.

#The computer is highly prone to error because the values it is computing are very small... i.e. the differences are tiny.
#The error compounds as the size of the Hilbert matrix increases.

#calculates coefficients for polynomial interpolation
def coef(x,y):
    a = [] #need array to store values.
    n = len(x)
    for i in range(n):
            a.append(y[i])
    for j in range(1,n):
        for i in reversed(range(j,n)):
            #print i,j
            #print ((a[i]-a[i-1])/(x[i]-x[i-j]))
            a[i]=((a[i]-a[i-1])/(x[i]-x[i-j])) #divided differences
    return a


def eval(x,y,p): #single real value. Function returns value of the interpolating polynomial at p.
    a = coef(x,y) #output form coef function.
    #print a #diagnostic check to see if a contains the correct values.
    n = len(x)
    #print n #diagnostic check to see if the length matches the number of values of x.
    t = a[n-1]
    for i in reversed(range(n-1)): #the loop that generates the polynomial.
        t = t*(p-x[i])+a[i]
    return t #returns the value.
  
  
#chebyshev x-value calculation. 

def chebyshevnodes(func,a,b,n): #number of nodes defined by user as n, a=left most point and b=right most point.
    x = [] #empty matrix to store values.
    y = []
    for i in range(1): #really only need to check this once.
        if 0<=i<=n: #diagnostic to make sure we are operating within the bounds of the function. Not really needed.
            print "i is less than or equal to n"
        else:
            print "i is NOT less than or equal to n"
            break
    #for i in range(n+1):
     #   x.append((0.5*(a+b))+((0.5*(b-a))*cos((2.0*float(i))/(2.0*float(n))*pi))) #formula for Chebyshev x-value calculation for ANY range [a,b].
    for i in range(n+1):
        x.append(cos((((2.0*float(i)+1)/(2.0*float(n)+2))*pi))) #formula given in book.
        y.append(func(x[i]))
    return [x,y] #returns calculated chebyshev nodes to matrix for use.

print chebyshevnodes(-5,5,10) #print the result. 

#Problem 4.1.1

xvalues = array([1,2,3,-4,5]) #x-values
yvalues = array([2,48,272,1182,2262]) #y-values

def coef(x,y):
    a = [] #need array to store values.
    n = len(x)
    for i in range(n):
            a.append(y[i])
    for j in range(1,n):
        for i in reversed(range(j,n)):
            #print i,j
            #print ((a[i]-a[i-1])/(x[i]-x[i-j]))
            a[i]=((a[i]-a[i-1])/(x[i]-x[i-j])) #divided differences
    return a

print coef(xvalues,yvalues)

def eval(x,y,p): #single real value. Function returns value of the interpolating polynomial at p.
    a = coef(x,y) #output form coef function.
    #print a #diagnostic check to see if a contains the correct values.
    n = len(x)
    #print n #diagnostic check to see if the length matches the number of values of x.
    t = a[n-1]
    for i in reversed(range(n-1)): #the loop that generates the polynomial.
        t = t*(p-x[i])+a[i]
    return t #returns the value.

print eval(xvalues,yvalues,-1) #returns 12.

#so coef and eval functions work. 

#Problem 4.2.1
#Find the interpolating polynomial of degree 20 for f(x)=1/((x^2)+1)

#need to first generate a table of 21 equally spaced nodes from [-5,5].

def f(x): #our function that will generate values for y for a given x.
    return 1/((x**2)+1)

def node_calc(func, a, b, n): #n=number of nodes desired minus 1, a=left most point, b=right most point, func=function.
    distance = abs(b-a) #need the distance between a and b.
    #print distance #check to see if distance from b to zero is correct
    spacing = float(distance)/n #need to use float for integer division.
    #print spacing #diagnostic check to see if we are assigning the correct spacing.
    #print "above value indicates the spacing between individual values of x."
    x = []
    y = []
    for i in range(n+1):
        x.append(a+i*spacing)
    for i in range(n+1):
        y.append(func(x[i]))
    return [x,y]


table =  node_calc(f,-5,5,20)
table2 = node_calc(f,-10,10,40)
#print table
#print table2

xvalues = table[0]
yvalues = table[1]
xvalues2 = table2[0]
yvalues2 = table2[1]

#print xvalues
#print yvalues
#print xvalues2
#print yvalues2

#need to define both coef and eval for print.

def coef(x,y):
    a = [] #need array to store values.
    n = len(x)
    for i in range(n):
            a.append(y[i])
    for j in range(1,n):
        for i in reversed(range(j,n)):
            #print i,j
            #print ((a[i]-a[i-1])/(x[i]-x[i-j]))
            a[i]=((a[i]-a[i-1])/(x[i]-x[i-j])) #divided differences
    return a

def eval(x,y,p): #single real value. Function returns value of the interpolating polynomial at p.
    a = coef(x,y) #output form coef function.
    #print a #diagnostic check to see if a contains the correct values.
    n = len(x)
    #print n #diagnostic check to see if the length matches the number of values of x.
    t = a[n-1]
    for i in reversed(range(n-1)): #the loop that generates the polynomial.
        t = t*(p-x[i])+a[i]
    return t #returns the value.

print "{:^50}\n".format("f(x) = 1/(x^2 + 1)")

print " {:^7} {:^14} {:^14} {:^14}\n".format('x', 'f(x)','p(x)','f(x)-p(x)')

for i in table[0]: #doing this to set a length for the loop and thus the correct amount to print.
    print " {:< 6.1f}".format(i),
    print " {:< 13e}".format(eval(xvalues2,yvalues2,i)), #evaluating f(x).
    print " {:< 13e}".format(f(i)), #evaluating p(x) using given x-values.
    print " {:< 13e}".format((f(i) - eval(xvalues2,yvalues2,i))) #error


#All over the output range, there is a huge discrepancy between f(x) and p(x). This shows just
#how easily polynomial interpolation can fail to approximate the polynomial at the endpoints, especially
#when a large amount of nodes are used. The values of the error are increasingly growing in magnitude as x increases,
#in a sort of wave pattern.

#Problem 4.2.2a

def f(x): #our function that will generate values for y for a given x.
    return 1/((x**2)+1)

def node_calc2(func,n): #n=number of nodes desired minus 1, a=left most point, b=right most point, func=function.
    x = []
    y = []
    for i in range(n+1):
        x.append(5*cos((i*pi)/n)) #chebyshev node calculation.
        y.append(func(x[i]))
    return [x,y]

table = node_calc2(f,40)
#print table

xvalues = table[0]
yvalues = table[1]

#for i in xvalues:
 #   print i

#print len(xvalues)

#print xvalues
#print yvalues2

#need to define both coef and eval for print.

def coef(x,y):
    a = [] #need array to store values.
    n = len(x)
    for i in range(n):
            a.append(y[i])
    for j in range(1,n):
        for i in reversed(range(j,n)):
            #print ((a[i]-a[i-1])/(x[i]-x[i-j]))
            a[i]=((a[i]-a[i-1])/(x[i]-x[i-j])) #divided differences
    return a


#print coef(xvalues2,yvalues2)


def eval(x,y,p): #single real value. Function returns value of the interpolating polynomial at p.
    a = coef(x,y) #output form coef function.
    #print a #diagnostic check to see if a contains the correct values.
    n = len(x)
    #print n #diagnostic check to see if the length matches the number of values of x.
    t = a[n-1]
    for i in reversed(range(n-1)): #the loop that generates the polynomial.
        t = t*(p-x[i])+a[i]
    return t #returns the value.

#print eval(xvalues,yvalues,-1)

print "{:^50}\n".format("f(x) = 1/(x^2 + 1)")

print " {:^7} {:^14} {:^14} {:^14}\n".format('x', 'f(x)','p(x)','f(x)-p(x)')

for i in table[0]: #doing this to set a length for the loop and thus the correct amount to print.
    print " {:< 6.1f}".format(i),
    print " {:< 13e}".format(eval(xvalues,yvalues,i)), #evaluating f(x).
    print " {:< 13e}".format(f(i)), #evaluating p(x) using given x-values.
    print " {:< 13e}".format((f(i) - eval(xvalues,yvalues,i))) #error


#Problem 4.2.2b

def f(x): #our function that will generate values for y for a given x.
    return 1/((x**2)+1)

def node_calc3(func,n): #n=number of nodes desired minus 1, a=left most point, b=right most point, func=function.
    x = []
    y = []
    for i in range(n+1):
        x.append(5*cos(((2*i+1.0)*pi)/((2*(n+1))))) #chebyshev node calculation.
        y.append(func(x[i]))
    return [x,y]

table = node_calc3(f,40)
#print table

xvalues = table[0]
yvalues = table[1]

#for i in xvalues:
 #   print i

#print len(xvalues)

print xvalues
#print yvalues

#need to define both coef and eval for print.

def coef(x,y):
    a = [] #need array to store values.
    n = len(x)
    for i in range(n):
            a.append(y[i])
    for j in range(1,n):
        for i in reversed(range(j,n)):
            #print i,j
            #print ((a[i]-a[i-1])/(x[i]-x[i-j]))
            a[i]=((a[i]-a[i-1])/(x[i]-x[i-j])) #divided differences
    return a


#print coef(xvalues,yvalues)

def eval(x,y,p): #single real value. Function returns value of the interpolating polynomial at p.
    a = coef(x,y) #output form coef function.
    #print a #diagnostic check to see if a contains the correct values.
    n = len(x)
    #print n #diagnostic check to see if the length matches the number of values of x.
    t = a[n-1]
    for i in reversed(range(n-1)): #the loop that generates the polynomial.
        t = t*(p-x[i])+a[i]
    return t #returns the value.

print "{:^50}\n".format("f(x) = 1/(x^2 + 1)")

print " {:^7} {:^14} {:^14} {:^14}\n".format('x', 'f(x)','p(x)','f(x)-p(x)')

for i in table[0]: #doing this to set a length for the loop and thus the correct amount to print.
    print " {:< 6.1f}".format(i),
    print " {:< 13e}".format(eval(xvalues,yvalues,i)), #evaluating f(x).
    print " {:< 13e}".format(f(i)), #evaluating p(x) using given x-values.
    print " {:< 13e}".format((f(i) - eval(xvalues,yvalues,i))) #error


#Problem 4.2.10

#chebyshev x-value calculation.

#need to define both coef and eval for print.

def coef(x,y):
    a = [] #need array to store values.
    n = len(x)
    for i in range(n):
            a.append(y[i])
    for j in range(1,n):
        for i in reversed(range(j,n)):
            a[i]=((a[i]-a[i-1])/(x[i]-x[i-j])) #divided differences
    return a

def eval(x,y,p): #single real value. Function returns value of the interpolating polynomial at p.
    a = coef(x,y) #output form coef function.
    #print a #diagnostic check to see if a contains the correct values.
    n = len(x)
    #print n #diagnostic check to see if the length matches the number of values of x.
    t = a[n-1]
    for i in reversed(range(n-1)): #the loop that generates the polynomial.
        t = t*(p-x[i])+a[i]
    return t #returns the value.

def f(x): #our function that will generate values for y for a given x.
    return 1/((x**2)+1)

def chebyshevnodes(func,a,b,n): #number of nodes defined by user as n, a=left most point and b=right most point.
    x = [] #empty matrix to store values.
    y = []
    for i in range(1): #really only need to check this once.
        if 0<=i<=n: #diagnostic to make sure we are operating within the bounds of the function. Not really needed.
            print "i is less than or equal to n"
        else:
            print "i is NOT less than or equal to n"
            break
    #for i in range(n+1):
     #   x.append((0.5*(a+b))+((0.5*(b-a))*cos((2.0*float(i))/(2.0*float(n))*pi))) #formula for Chebyshev x-value calculation for ANY range [a,b].
    for i in range(n+1):
        x.append(cos((((2.0*float(i)+1)/(2.0*float(n)+2))*pi))) #formula given in book.
        y.append(func(x[i]))
    return [x,y] #returns calculated chebyshev nodes to matrix for use.

#chvvalues41 = chebyshevnodes(-1,1,40) #needs to be 40 because loop goes to n+1.
#print chvvalues41

#print "{:^50}\n".format("f(x) = 1/(x^2 + 1)")

#print " {:^7} {:^14} {:^14} {:^14}\n".format('x', 'f(x)','p(x)','f(x)-p(x)')

#for i in table[0]: #doing this to set a length for the loop and thus the correct amount to print.
 #   print " {:< 6.1f}".format(i),
  #  print " {:< 13e}".format(eval(chvvalues41,yvalues,i)),
   # print " {:< 13e}".format(f(i)),
    #print " {:< 13e}".format((f(i) - eval(chvvalues41,yvalues,i)))

#print chebyshevnodes(-1,1,3) #print the result.
#print chebyshevnodes(-1,1,7) #print the result.
#print chebyshevnodes(-1,1,15) #print the result.

chvvalues15 = chebyshevnodes(f,-1,1,15)

xvalues = chvvalues15[0]
yvalues = chvvalues15[1]

#print xvalues

a = -1
b = 1

cap = 2.**-15.
#print cap

print "\n {:^5.1}".format('x'), " {:^12}".format('p(x)'), " {:^12}".format('2^-n'), " {:^12}".format('2^n-1 - p(x)')

for i in chvvalues15[0]:
    print " {:< 5.1f}".format(i),
    print " {:<12.6e}".format(eval(xvalues,yvalues,i)),
    print " {:<12.6e}".format(cap),
    print " {:< 12.6e}".format((cap) - eval(xvalues,yvalues,i))


#chebyshev nodes are better because it reduces the oscillations normally associated with polynomial interpolation.
#We can clearly see here that the error never exceeds 2.**-15.

def deriv(f,x):
    n = 5
    h = 1.
    d = []
    for i in range(n+1):
        d.append([(f(x + h) - f(x - h))/(2. * h)])
        for j in range(1,i+1):
            d[i].append(d[i][j-1]+(d[i][j-1]-d[i-1][j-1])/(4.0**j - 1.0))
        h /= 2.
    return d
    
def deriv2(f,x):
    n = 5
    h = 1.
    d = []
    for i in range(n+1):
        d.append([(f(x + h) - 2. * f(x) + f(x - h))/(h ** 2.)])
        for j in range(1,i+1):
            d[i].append(d[i][j-1]+(d[i][j-1]-d[i-1][j-1])/(4.**(j+1) - 1.))
        h /= 2.
    return d[n][n]

#Problem 4.3.1

def deriv(f,x):
    n = 5
    h = 1.
    d = []
    for i in range(n+1):
        #print i
        #print (f(x + h) - f(x - h))/(2. * h)
        d.append([(f(x + h) - f(x - h))/(2. * h)])
        for j in range(1,i+1):
            d[i].append(d[i][j-1]+(d[i][j-1]-d[i-1][j-1])/(4.0**j - 1.0))
        h /= 2.
    return d[n][n]


cos = lambda x: math.cos(x)
arctan = lambda x: math.atan(x)
abss = lambda x: abs(x)
#print deriv(cos, 0.)
#print deriv(arctan, 1.)
#print deriv(abss, 0.)

def resultsmatrix(n,f,x): #function to create matrix A
    A = np.empty((n,n), order='F') #declare an empty array of nxn dimension
    for i in range(n):
        for j in range(n):
            A[i,j]=deriv(f,x)
    print "The value of the derivative at x is {:.5}".format(A[1,1]) #print derivative from element 1,1.
    return A

print resultsmatrix(10,cos,0.)

print resultsmatrix(10,arctan,1.)

print resultsmatrix(10,abss,0.)

#Problem 4.3.2

def deriv2(f,x):
    n = 5
    h = 1.
    d = []
    for i in range(n+1):
        d.append([(f(x + h) - 2. * f(x) + f(x - h))/(h ** 2.)])
        for j in range(1,i+1):
            d[i].append(d[i][j-1]+(d[i][j-1]-d[i-1][j-1])/(4.**(j+1) - 1.))
        h /= 2.
    return d[n][n]


cos = lambda x: math.cos(x)
arctan = lambda x: math.atan(x)
abss = lambda x: abs(x)

def resultsmatrix2(n,f,x): #function to create matrix A
    A = np.empty((n,n), order='F') #declare an empty array of nxn dimension
    for i in range(n):
        for j in range(n):
            A[i,j]=deriv2(f,x)
    print "The value of the 2nd derivative at x is {:.5}".format(A[1,1]) #print derivative from element 1,1.
    return A

print resultsmatrix2(10,cos,0.)

print resultsmatrix2(10,arctan,1.)

print resultsmatrix2(10,abss,0.)


#____________#



def deriv(f, x, h, n):
    d = np.matrix(np.zeros((n, n)))
    for i in range(0, n):
        d[i, 0] = (f(x + h) - f(x - h)) / (2. * h)
        for j in range(1, i + 1):
            d[i, j] = d[i, j - 1] + (d[i, j - 1] - d[i - 1, j - 1]) / (4. ** j - 1.)
        h = h / 2.
    return d


my_cos = lambda x: math.cos(x)
my_arctan = lambda x: math.atan(x)
my_abs = lambda x: abs(x)
print deriv(my_cos, 0, 1., 5)
print deriv(my_arctan, 1., 1., 5)
print deriv(my_abs, 0, 1., 5)


def gen_coef(xlist,ylist):    #this function generates the "divided difference" coefficients (ai's)
    alist = []    #creating a list to store the coefficients
    n = len(xlist)    #setting the step number
    for i in range(n):    #starting with y values
        alist.append(ylist[i])
    for j in range (1,n):
        for i in reversed(range(j,n)):
            alist[i] = ((alist[i] - alist[i-1])/(xlist[i]-xlist[i-j]))    #calculating the successive divided differences
    return alist    #return the value


def inter_eval(xlist,ylist,value):    #this evaluates values using the interpolating polynomial for a set of x's and y's
    alist = gen_coef(xlist,ylist)    #uses gen_coef to generate coefficients
    n = len(xlist)    #setting the step number
    temp = alist[n-1]    #starting with the last coefficient for a horners algorithm style evalutation
    for i in reversed(range(n-1)):    #evaluation
        temp *= (value - xlist[i])
        temp += alist[i]
    return temp    #returning the value

xlist = [1., 2., 3., -4., 5.]
ylist = [2., 48., 272., 1182., 2262.]

print "\nFrom the table:\n"
print "x",
for i in range(len(xlist)): print "| {:<7}".format(xlist[i]),
print
print "y",
for i in range(len(ylist)): print "| {:<7}".format(ylist[i]),
print
print "\nThe interpolating polynomial p(x) was determined.\n"
print "{:^50}".format("p(-1) = " + str(inter_eval(xlist,ylist, -1.)))

#problem 1 on 4.2

def gen_table(func, start, stop, steps):
    step_size = (stop - start)/steps
    xlist = []
    ylist = []
    for i in range(steps + 1):
        xlist.append(start + i * step_size)
    for i in range(steps + 1):
        ylist.append(f(xlist[i]))
    return (xlist,ylist)
    

def gen_table2(func, steps):
    xlist = []
    ylist = []
    for i in range(steps + 1):
        xlist.append(5*math.cos((i*math.pi)/20.0))
    for i in range(steps + 1):
        ylist.append(f(xlist[i]))
    return (xlist,ylist)
    
def gen_table3(func, steps):
    xlist = []
    ylist = []
    for i in range(steps + 1):
        xlist.append(5*math.cos(((2*i+1.0)*math.pi)/42.0))
    for i in range(steps + 1):
        ylist.append(f(xlist[i]))
    return (xlist,ylist)


f = lambda x: 1.0/(x**2.0 + 1.0)


tables1 = gen_table(f, -5.0, 5.0, 20)

tables2 = gen_table(f, -10., 10., 40)

print """\nThe interpolating polynomial p(x) is created from a table that is
made by evaluating f(x) at 21 nodes of x that are evenly spaced
on the interval -5 <= x <= 5.

We then compared the difference between f(x) and p(x) at 41 x that were
evenly spaced on the interval -10.0 <= x <= 10.0.\n"""

print "{:^50}\n".format("f(x) = 1/(x^2 + 1)")

print " {:^7} {:^14} {:^14} {:^14}\n".format('x', 'f(x)','p(x)','f(x)-p(x)')

for i in tables2[0]:
    print " {:< 6.1f}".format(i),
    print " {:< 13e}".format(inter_eval(tables1[0],tables1[1],i)),
    print " {:< 13e}".format(f(i)),
    print " {:< 13e}".format((f(i) - inter_eval(tables1[0],tables1[1],i)))
    
tables1 = gen_table2(f, 20)

print """\nThe preceding test was repeated, this time with the interpolating
polynomial p(x) computed using chebyshev nodes for the x values,
calculated with as follows:

            xi = 5cos((i*pi)/20)\n"""

print " {:^7} {:^14} {:^14} {:^14}\n".format('x', 'f(x)','p(x)','f(x)-p(x)')

for i in tables2[0]:
    print " {:< 6.1f}".format(i),
    print " {:< 13e}".format(inter_eval(tables1[0],tables1[1],i)),
    print " {:< 13e}".format(f(i)),
    print " {:< 13e}".format((f(i) - inter_eval(tables1[0],tables1[1],i)))

tables1 = gen_table3(f, 20)

print """\nThe preceding test was repeated, this time with the interpolating
polynomial p(x) computed using chebyshev nodes for the x values,
calculated with as follows:

            xi = 5cos((2i+1)*pi)/42)\n"""

print " {:^7} {:^14} {:^14} {:^14}\n".format('x', 'f(x)','p(x)','f(x)-p(x)')

for i in tables2[0]:
    print " {:< 6.1f}".format(i),
    print " {:< 13e}".format(inter_eval(tables1[0],tables1[1],i)),
    print " {:< 13e}".format(f(i)),
    print " {:< 13e}".format((f(i) - inter_eval(tables1[0],tables1[1],i)))    
    
#hmwk 4.2.10

def test_cheby(num, start, stop):
    xlist = []
    for i in range(num + 1):
        xlist.append(math.cos(((2.0*float(i)+1.0)*math.pi)/(2.0*float(num) + 2.0)))
    output = []
    for i in range(start,stop +1):
        temp = (float(i)/10.0 - xlist[0])   #starting with the last coefficient for a horners algorithm style evalutation
        for i in range(1,len(xlist)):    #evaluation
            temp *= (float(i)/10.0 - xlist[i])
        output.append(abs(temp)) 
    return output
    
    
start, stop = -10, 10    

print "\nFor n = 3, we get the following results from testing error values of chebyshev nodes\n"

cheby = test_cheby(3, start, stop)    

cap = 2.**-3.

print "\n {:^5.1}".format('x'), " {:^12}".format('p(x)'), " {:^12}".format('2^-n'), "  {:^12}".format('2^n-1 - p(x)')
    
for i in range(stop - start + 1):
    print " {:< 5.1f}".format(i+start), " {:<12.6e}".format(cheby[i]), " {:<12.6e}".format(cap), " {:< 12.6e}".format(cap - cheby[i])
 
print "\n\nFor n = 7, we get the following results from testing error values of chebyshev nodes\n"
   
cheby = test_cheby(7, start, stop)    

cap = 2.**-7.

print "\n {:^5.1}".format('x'), " {:^12}".format('p(x)'), " {:^12}".format('2^-n'), "  {:^12}".format('2^n-1 - p(x)')
    
for i in range(stop - start + 1):
    print " {:< 5.1f}".format(i+start), " {:<12.6e}".format(cheby[i]), " {:<12.6e}".format(cap), " {:< 12.6e}".format(cap - cheby[i])

print "\n\nFor n = 15, we get the following results from testing error values of chebyshev nodes\n"

cheby = test_cheby(15, start, stop)    

cap = 2.**-15.

print "\n {:^5.1}".format('x'), " {:^12}".format('p(x)'), " {:^12}".format('2^-n'), "  {:^12}".format('2^n-1 - p(x)')
    
for i in range(stop - start + 1):
    print " {:< 5.1f}".format(i+start), " {:<12.6e}".format(cheby[i]), " {:<12.6e}".format(cap), " {:< 12.6e}".format(cap - cheby[i])


#4.3.1

def deriv(f,x):
    n = 5
    h = 1.
    d = []
    for i in range(n+1):
        d.append([(f(x + h) - f(x - h))/(2. * h)])
        for j in range(1,i+1):
            d[i].append(d[i][j-1]+(d[i][j-1]-d[i-1][j-1])/(4.0**j - 1.0))
        h /= 2.
    return d[n][n]
    
def dblderiv(f,x):
    n = 5
    h = 1.
    d = []
    for i in range(n+1):
        d.append([(f(x + h) - 2. * f(x) + f(x - h))/(h ** 2.)])
        for j in range(1,i+1):
            d[i].append(d[i][j-1]+(d[i][j-1]-d[i-1][j-1])/(4.**(j+1) - 1.))
        h /= 2.
    return d[n][n]

my_cos = lambda x: math.cos(x)
my_arctan = lambda x: math.atan(x)
my_abs = lambda x: abs(x)
print
print "The derivitive of cos(x) at x = 0 is {:.5}".format(deriv(my_cos, 0))
print "The derivitive of arctan(x) at x = 1 is {:.5}".format(deriv(my_arctan, 1.))
print "The derivitive of abs(x) at x = 0 is {:.5}".format(deriv(my_abs, 0))
print

#4.3.2

print "The 2nd derivitive of cos(x) at x = 0 is {:.5}".format(dblderiv(my_cos, 0.))
print "The 2nd derivitive of cos(x) at x = 1 is {:.5}".format(dblderiv(my_arctan, 1.))
print "The 2nd derivitive of cos(x) at x = 0 is {:.5}".format(dblderiv(my_abs, 0.)) #error
def naive_gauss(a, b):
    n = a.shape[0] - 1
    x = [0 for i in range(n)]
    for k in range(0, n):
        print a, '\n'
        for i in range(k + 1, n + 1):
            xmult = a[i, k] / a[k, k]
            a[i, k] = xmult
            for j in range(k + 1, n + 1):
                if i == 1 and j == 0: print "i=1 j=0"
                a[i, j] = a[i, j] - xmult * a[k, j]
            b[i] = b[i] - xmult * b[k]

    x.insert(n, (b[n] / a[n, n]).item(0))

    print 'Now back substitute'
    print a
    print b
    for i in range(n, -1, -1):
        s = b[i]
        print
        for j in range(n, i, -1):
            s = s - a[i, j] * x[j]
        x[i] = s / a[i, i]
    return x


a = np.matrix('6. -2. 2. 4.;12. -8. 6. 10.;3. -13. 9. 3.;-6. 4. 1. -18.')
b = np.matrix('16. 26. -19. -34.').transpose()
print naive_gauss(a, b)

n = 4
a = np.matrix(np.zeros([n, n]))
b = np.matrix(np.zeros(n)).transpose()
for i in range(0, n):
    b[i] = i + 1
    for j in range(0, n):
        a[i, j] = i + j

print 'a\n', a
print '\nb\n', b

print naive_gauss(a, b)
print 'We would expect infinite solutions since we have the following equations:'
print 'x1   -x3 -2x4 = 0'
print 'x2 + 2x3 +3x4= 1'


def naive_complex_gauss(a, b):
    n = a.shape[0] - 1
    x = [0 for i in range(n)]
    for k in range(0, n):
        imag = -a[k, k].item().imag
        for l in range(k, n):
            if k == l:
                a[k, l] = a[k, k].item().real
            else:
                a[k, l] = a[k, l] + np.complex(0, imag)
                #print a
        for i in range(k + 1, n + 1):
            xmult = a[i, k] / a[k, k]
            a[i, k] = xmult
            for j in range(k + 1, n + 1):
                if i == 1 and j == 0: print "i=1 j=0"
                a[i, j] = a[i, j] - xmult * a[k, j]
            b[i] = b[i] - xmult * b[k]

    x.insert(n, (b[n] / a[n, n]).item(0))
    print a
    print 'Now back substitute'
    print
    for i in range(n, -1, -1):
        s = b[i]
        for j in range(n, i, -1):
            #print s
            s = s - a[i, j] * x[j]
        x[i] = s / a[i, i]
    return x


a = np.matrix([[np.complex(5, 9), np.complex(5, 5), np.complex(-6, -6), np.complex(-7, -7)],
               [np.complex(3, 3), np.complex(6, 10), np.complex(-5, -5), np.complex(-6, -6)],
               [np.complex(2, 2), np.complex(3, 3), np.complex(-1, 3), np.complex(-5, -5)],
               [np.complex(1, 1), np.complex(2, 2), np.complex(-3, -3), np.complex(0, 4)]])
b = np.matrix([np.complex(-10, 2), np.complex(-5, 1), np.complex(-5, 1), np.complex(-5, 1)]).transpose()
print naive_complex_gauss(a, b)
b = np.matrix([np.complex(2, 6), np.complex(4, 12), np.complex(2, 6), np.complex(2, 6)]).transpose()
print naive_complex_gauss(a, b)
b = np.matrix([np.complex(7, -3), np.complex(7, -3), np.complex(0, 0), np.complex(7, -3)]).transpose()
print naive_complex_gauss(a, b)
b = np.matrix([np.complex(-4, -8), np.complex(-4, -8), np.complex(-4, -8), np.complex(0, 0)]).transpose()
print naive_complex_gauss(a, b)


def naive_gauss(a, b):
    n = a.shape[0] - 1
    x = [0 for i in range(n)]
    for k in range(0, n):
        print a, '\n'
        for i in range(k + 1, n + 1):
            xmult = a[i, k] / a[k, k]
            a[i, k] = xmult
            for j in range(k + 1, n + 1):
                if i == 1 and j == 0: print "i=1 j=0"
                a[i, j] = a[i, j] - xmult * a[k, j]
            b[i] = b[i] - xmult * b[k]
    x.insert(n, (b[n] / a[n, n]).item(0))
    print 'Now back substitute'
    print a
    print b
    for i in range(n, -1, -1):
        s = b[i]
        print
        for j in range(n, i, -1):
            s = s - a[i, j] * x[j]
        x[i] = s / a[i, i]
    return x


def gauss(n, a):
    s = [0 for i in range(n)]
    x = [0 for i in range(n)]
    l = [i for i in range(n)]
    for i in range(0, n):
        l[i] = i
        smax = 0
        for j in range(0, n):
            #print abs(a[i,j])
            smax = max(smax, abs(a[i, j]))
        s[i] = smax
        #print s
        #print
    for k in range(0, n - 1):
        rmax = 0
        for i in range(k, n):
            #print k/s[l[i]]
            r = abs(a[l[i], l[k]] / s[l[i]])
            #print r>rmax
            #print 'j=',j,' i=',i, 'r=',r,'rmax=',rmax
            if r > rmax:
                rmax = r
                j = i

        #print k
        temp = l[k]
        l[k] = l[j]
        l[j] = temp
        #print l
        #print '\n'
        for i in range(n - 1, k, -1):
            xmult = a[l[i], k] / a[l[k], k]
            #print xmult
            a[l[i], k] = xmult
            for j in range(k + 1, n):
                a[l[i], j] = a[l[i], j] - xmult * a[l[k], j]
    return a, l


a = np.matrix('.4096 .1234 .3678 .2943;' +
              '.2246 .3872 .4015 .1129;' +
              '.3645 .1920 .3781 .0643;' +
              '.1784 .4002 .2786 .3927')
n = 4
a, l = gauss(n, a)
print
print a
print
print l
print

b = [.4043, .1550, .4240, .2557]


def solve(n, a, l, b):
    x = [0 for i in range(n)]
    for k in range(0, n - 1):
        for i in range(k + 1, n):
            b[l[i]] = b[l[i]] - a[l[i], k] * b[l[k]]
            #print n,x
    x[n - 1] = b[l[n - 1]] / a[l[n - 1], n - 1]
    for i in range(n - 1, -1, -1):
        s = b[l[i]]
        for j in range(i + 1, n):
            s = s - a[l[i], j] * x[j]
        x[i] = s / a[l[i], i]
    return x


print solve(n, a, l, b)


def find_b(n):
    b = [0 for i in range(n)]
    for i in range(0, n):
        b[i] = n * (i + 1.) * ((i + 1.) + n + 1.) + (1. / 6) * n * (1. + n * (2. * n + 3.))
    return b


def find_a(n):
    a = np.matrix(np.zeros([n, n]))
    for i in range(1, n + 1):
        for j in range(1, n + 1):
            a[i - 1, j - 1] = (i + j) ** 2
    return a


n = 2
a1 = find_a(n)
b1 = find_b(n)
print a1
print b1
a1, l = gauss(n, a1)
print solve(n, a1, l, b1)

print '\nn=3'
n = 3
a2 = find_a(n)
b2 = find_b(n)
print a2
print b2
a2, l = gauss(n, a2)
print solve(n, a2, l, b2)

print '\nn=4'
n = 4
a3 = find_a(n)
b3 = find_b(n)
print a3
print b3
a3, l = gauss(n, a3)
print solve(n, a3, l, b3)

a = np.matrix('.0001 -5.03  5.809 7.832;' +
              '2.266 1.995  1.212 8.008;' +
              '8.850 5.681  4.552 1.302;' +
              '6.775 -2.253 2.908 3.970')

b = np.matrix('9.574 7.219 5.73 6.291').transpose()
n = 4
print '\nUsing naive_gauss'
print naive_gauss(a, b)
a, l = gauss(n, a)
print '\nUsing gauss and solve'
print solve(n, a, l, b)

print """
def trapezoid_uniform(f,a,b,n):
    h=(b-a)/n
    s=.5*(f(a)+f(b))
    for i in range(1,n-1):
        x=a+i*h
        s=s+f(x)
    return s*h
"""


def trapezoid_uniform(f, a, b, n):
    h = (b - a) / n
    s = .5 * (f(a) + f(b))
    for i in range(1, n - 1):
        x = a + i * h
        s = s + f(x)
    return s * h


def f_1(x): return math.sin(x)


def f_2(x): return math.e ** x


def f_3(x): return math.atan(x)


print trapezoid_uniform(f_1, 0, math.pi, 100)
print trapezoid_uniform(f_2, 0, 1., 100)
print trapezoid_uniform(f_3, 0, 1., 100)

low_high = [.1, 1., 10., 100., 10000., 1e60]


def calc_low_highs(f, lh):
    for b in lh:
        print 'a = 0 b=', b
        try:
            print trapezoid_uniform(f, 0.5, b, 100)
            print
        except OverflowError:
            print 'Overflow error, unable to compute\n'


def f_4(x): return 1. / (math.e ** (x ** 2))


def f_5(x): return (1. / x) * math.sin(x)


def f_6(x): return math.sin(x ** 2)


calc_low_highs(f_4, low_high)
calc_low_highs(f_6, low_high)


def f_7(x): return 1. / math.e ** (-math.log(x ** 2))


def f_8(x): return x * math.sin(x)


def f_9(x): return math.sin(math.tan(x ** 2))


calc_low_highs(f_7, low_high)
calc_low_highs(f_8, low_high)
calc_low_highs(f_9, low_high)

#calculates f(x) on the interval [a,b] using the trapezoid rule with n equal sub-intervals.
def trapezoid_uniform(f,a,b,n): #f denotes f(x), a=left most point, b=right most point, n=number of sub-intervals.
    h = (b-a)/n #sub-interval size
    sum = 0.5*(f(a)+f(b))
    for i in range(n-1):
        x = a+float(i)*h
        sum = sum + f(x)
    sum = sum*h
    return sum #output

#Test functions

def function1(x):
    return sin(x)

def function2(x):
    return e**x

def function3(x):
    return atan(x)
   
    
print trapezoid_uniform(function1,0.,pi,100) #answer is 2, we have 1.99884870579.
print trapezoid_uniform(function2,0.,1.,100) #answer is 1.7828, we have 1.70138380273.
print trapezoid_uniform(function3,0.,1.,100) #answer is 0.438825, we have 0.431016675615.








def simpson(f, a, b, level, level_max, p=1e-5):
    level += 1
    h = b - a
    c = (a + b) / 2
    one_s = h * (f(a) + 4. * f(c) + f(b)) / 6.
    d = (a + c) / 2
    e = (c + b) / 2
    two_s = h * (f(a) + 4. * f(d) + 2. * f(c) + 4. * f(e) + f(b)) / 12.
    if level >= level_max:
        simpson_result = two_s
        print 'Max level reached'
        return simpson_result
    else:
        if abs(two_s - one_s) < 15. * p:
            return two_s + (two_s - one_s) / 15
        else:
            left_s = simpson(f, a, c, level, level_max, p / 2.)
            right_s = simpson(f, c, b, level, level_max, p / 2.)
            print left_s
            return left_s + right_s


def f_1(x): return math.cos(2. * x) / math.e ** x


def f_2(x): return 1 / (1 + x ** 2)


def f_3(x): return math.sqrt(1 - x ** 2) - x


def f_4(x): return math.cos(2. * x) / math.e ** x


a = 0
b = 1.

print 4 * simpson(f_2, a, b, 0, 5)

a = 0
b = 1. / math.sqrt(2)
print 8 * simpson(f_3, a, b, 0, 5)

a = 0
b = 5. * math.pi / 4.
print simpson(f_4, a, b, 0, 5)
