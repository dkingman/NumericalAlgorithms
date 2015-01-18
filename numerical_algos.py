# -*- coding: utf-8 -*-
"""
Created on Tue Jan 29 23:14:17 2013

@authors: David Kingman, Valentino Constantinou, Taylor Stevens
"""

#NOTE: There are currently some errors in the code (linear algebra). 

import math
import copy
import matplotlib.pyplot as plt
from numarray import argmax
from math import cos
from math import tan
from math import atan
from math import sin
from math import pi
from math import e
from math import sqrt
from numpy import array
from numpy import zeros
from numpy import linalg

#distance between two points. Can easily extend to R^n. 
def point_distance_R2(x1, y1, x2, y2):
    return sqrt(((x2 - x1) ** 2.) + ((y2 - y1) ** 2.))

def point_distance_R3(x1, y1, z1, x2, y2, z2):
    return sqrt(((x2 - x1) ** 2.) + ((y2 - y1) ** 2.) + ((z2 - z1) ** 2.))

print point_distance_R2(0., 1., 1., 0.) #Pythagorean Theorem tells us this is the correct result.
print point_distance_R3(1., 2., 3., 3., 2., 1.) #correct according to Wolfram Alpha

#calculate the hypotenuse of a triangle
def hypot(x, y):
    x = abs(x)
    y = abs(y)
    if x > y:
        r = y / x
        return x * sqrt(1 + r * r)
    if y == 0:
        return 0.
    r = x / y
    return y * sqrt(1 + r * r)

print hypot(3., 4.) #returns 5

#area of a triangle using points
def tri_point_area(x1,x2,x3,y1,y2,y3):
    A = [[x1,y1,1.],[x2,y2,1.],[x3,y3,1.]] #put vector of x's, vector of y's, and vector of 1's as columns in 3x3 matrix.
    #print A
    area = 0.5*abs(linalg.det(A))
    return area

print tri_point_area(-1.,0.,1.,0.,1.,0.)


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

def newton(f, f_prime, x, tolerance, precision, display, steps):
    fx= f(x)
    for i in range(steps):
        if abs(f_prime(x)) < precision:
            print 'small derivative'
            break
        d = fx / f_prime(x)
        x = x - d
        fx = f(x)
        if display == 1:
                print [i, x, fx]
        if abs(d) < tolerance:
            print 'convergence'
            break


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
    fx = f(x)
    for i in range(steps):
        if abs(f_prime(x)) < precision:
            print 'small derivative'
            break
        d = fx / f_prime(x)
        x = x - d
        if abs(f(x - d)) >= abs(f(x)):
            d = 0.5 * d
        else:
            d = d
        fx = f(x)
        if display == 1:
                print [i, x, fx]
        #print [steps,x,fx]
        if abs(d) < tolerance:
            print 'convergence'
            break



def testfn(x):
    return sin(x)

def testfnprime(x):
    return cos(x)


print newton(testfn,testfnprime,1.2,1*e**-6,1*e**-3,1,25)
print newtonmod(testfn,testfnprime,1.2,1*e**-6,1*e**-3,1,25)

#yay it works.

def newtonaccel(f, f_prime, x, tolerance, precision, display, steps):
    fx= f(x)
    for i in range(steps):
        if abs(f_prime(x)) < precision:
            print 'small derivative'
            break
        d = fx / f_prime(x)
        x = x - (2*d)
        fx = f(x)
        if display == 1:
                print [i, x, fx]
        if abs(d) < tolerance:
            print 'convergence'
            break


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
            print [n, a, fa]
        a = a - d
        fa = f(a)


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

#######
def createA(n): #function to create matrix A
    A = zeros((n,n), order='F') #declare an empty array of nxn dimension
    for i in range(n):
        for j in range(n):
             if j>=1:
                 A[1,j]=A[j,1]=n**(-1)
             else:
                 A[i,j]=A[i-1,j]+A[i,j-1]
    return A

A = createA(25)
#print A

def transpose_matrix(A): #where A is a matrix
    AT = []
    for i in range(len(A[0])):
        line = []
        for j in range(len(A)):
            line.append(A[j][i])
        AT.append(line)
    return AT #returns a transposed matrix.

A = [[1,2,3],[4,5,6]]
print transpose_matrix(A)


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
    A = zeros((n,n), order='F') #declare an empty array of nxn dimension
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
    A = zeros((n,n), order='F') #declare an empty array of nxn dimension
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
    A = zeros((n,n), order='F') #declare an empty array of nxn dimension
    for i in range(n):
        for j in range(n):
             A[i,j]= -1+2*argmax([i,j])
    return A

A = createA(30)
print A

def createB(n): #function to create matrix B
    B = zeros((n,1), order='F') #declare an empty array of nxn dimension
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
    H = zeros((n, n), order='F') #declare an empty array of nxn dimension.
    for j in range(n):
        for i in range(n):
            H[i,j] = 1. / (i + j - 1) #formula for constructing matrix values.
    return H #return the matrix A

H3 = createHilbert(3) #create 3x3 Hilbert matrix
print H3

H8 = createHilbert(8) #create 8x8 Hilbert matrix
#print H8

H13 = createHilbert(13) #create 13x13 Hilbert matrix
#print H13

def createB(n): #function to create matrix b
    H = zeros((n, n), order='F') #declare an empty array of nxn dimension.
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

#Cramer'srule for finding coefficients of linear systems.

#We have the system Ax=b, where A is an nxn matrix and b is an 1xn matrix. There are several ways to solve the system,
#but the algorithm below evaluates x1, x2, ..., xn using Cramer's Rule.

def cramers_rule_R2(a1, a2, a3, a4, b1, b2):
    A = [[a1, a2], [a3, a4]]
    if linalg.det(A) == 0:
        print "Cramer's rule does not apply."
    else:
        A1 = [[b1, a2], [b2, a4]]  #replace the first column with the matrix b.
        A2 = [[a1, b1], [a3, b2]]  #replace the second column with the matrix b.
        x1 = linalg.det(A1) / linalg.det(A)
        x2 = linalg.det(A2) / linalg.det(A)
        return [x1, x2]

def cramers_rule_R3(a1, a2, a3, a4, a5, a6, a7, a8, a9, b1, b2, b3):
    A = [[a1, a2, a3], [a4, a5, a6], [a7, a8, a9]]
    if linalg.det(A) == 0:
        print "Cramer's rule does not apply."
    else:
        A1 = [[b1, a2, a3], [b2, a5, a6], [b3, a8, a9]]  #replace the first column with the matrix b.
        A2 = [[a1, b1, a3], [a4, b2, a6], [a7, b3, a9]]  #replace the second column with the matrix b.
        A3 = [[a1, a2, b1], [a4, a5, b2], [a7, a8, b3]]  #replace the third column with the matrix b.
        x1 = linalg.det(A1) / linalg.det(A)
        x2 = linalg.det(A2) / linalg.det(A)
        x3 = linalg.det(A3) / linalg.det(A)
        return [x1, x2, x3]

print cramers_rule_R2(7., -2., 3., 1., 3., 5.)
print cramers_rule_R3(1., 0., 2., -3., 4., 6., -1., -2., 3., 6., 30., 8.)


#Polynomial Interpolation

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
  
#Node Transformations.
#Any arbitrary finite interval can be transformed to the interval [-1,1] using the first function provided below.

#transforming points from [a,b] to [-1,1]
def transform_1(array,a,b):
    x = [] #empty matrix to store values
    for i in range(len(array)):
        x.append(((2. * array[i]) - b - a) / (b - a))
    return x

#This next function is simply the inverse of the first transform function, which takes nodes in an interval [a.b] and
#transforms them to any arbitrary finite interval [a,b].

#transforming points from [-1,1] to [a,b]
def transform_2(array,a,b):
    x = [] #empty matrix to store values
    for i in range(len(array)):
        x.append(((array[i] * (b - a)) + (b + a)) / 2.)
    return x

z = [-4., -3., -2., -1., 0., 1., 2., 3., 4.]
z2 = transform_1(z, -4, 4)
print z2

print transform_2(z2, -4, 4)

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
    for i in range(1,n):
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
   
    
print trapezoid_uniform(function1,0.,pi,20000) #answer is 2, we have 1.99999999589
print trapezoid_uniform(function2,0.,1.,20000) #answer is 1.7828, we have 1.71828182882
print trapezoid_uniform(function3,0.,1.,20000) #answer is 0.438825, we have 0.438824573013

#trapazoid method, while generally affective, does not estimate the area under f(x)=e^x well. This is due to the
#exponentially increasing nature of the function. The higher the magnitude of the slope, the more error we can
#expect from using the trapezoid method.

#Gaussian probability integral
def f1(x):
    return e**-(x**2)

def f1b(t):
    return (1/t)*e**-(log(t))**2

#sine integral
def f2(x):
    return sin(x)/x

def f2b(t):
    return sin(1/t)/t

#Fresnel sine integral
def f3(x):
    return sin(x**2)

def f3b(t):
    return sin(tan(t)**2)/(cos(t)**2)

print trapezoid_uniform(f1,0.,200.,20000) #0.896226925453
print trapezoid_uniform(f1b,(10.**-10.),1.,20000) #returns 0.886176922744, which is close to the x-form.

print trapezoid_uniform(f2,(10.**-10.),200.,20000) #returns 1.56842626819, the correct result is 1.57
print trapezoid_uniform(f2b,(10.**-10.),200.,20000) #returns -24375299.9242

print trapezoid_uniform(f3,0.,200.,20000) #returns 0.631141945821, which is close to the true result of 0.626657
print trapezoid_uniform(f3b,0.,((0.5*pi)-(10.**-10.)),20000) #returns -3.84967571979e+15.


#______________#

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

#Open formulas are slightly better when only two or three points are used. When using more than three points,
#closed formulas are far more accurate than the open formulas. Additionally, a high-order formula may produce larger
#error than a low-order one. As a general rule, formulas employing more than eight points are almost never used.

#open Newton-Cotes Rules
def nc_midpoint(f,a,b,n):
    h = 0.5*(b-a) #midpoint
    x1 = a+h #we need to move one step away from a. 
    sum = 2.*h*f(x1) #formula
    return sum

def nc_two_point(f,a,b,n):
    h = 0.3333333*(b-a) #divide interval by 3
    x1 = a+h
    x2 = a+2.*h #we need to move twice the step size away from a.
    sum = 1.5*h*(f(x1)+f(x2)) #formula
    return sum

def nc_three_point(f,a,b,n):
    h = 0.25*(b-a) #divide interval by 4
    x1 = a+h #we need to move one step away from a. 
    x2 = a+2.*h #we need to move twice the step size away from a.
    x3 = a+3.*h #we need to move three times the step size away from a.
    sum = 1.33333333*h*(2.*f(x1)-f(x2)+2.*f(x3)) #formula
    return sum

def nc_four_point(f,a,b,n):
    h = 0.2*(b-a) #divide interval by 5
    x1 = a+h #we need to move one step away from a. 
    x2 = a+2.*h #we need to move twice the step size away from a.
    x3 = a+3.*h #we need to move three times the step size away from a.
    x4 = a+4.*h #we need to move four times the step size away from a.
    sum = 0.208333333*h*(11.*f(x1)+f(x2)+f(x3)+11.*f(x4)) #formula
    return sum

def nc_five_point(f,a,b,n):
    h = 0.1666666667*(b-a) #divide interval by 6
    x1 = a+h #we need to move one step away from a. 
    x2 = a+2.*h #we need to move twice the step size away from a.
    x3 = a+3.*h #we need to move three times the step size away from a.
    x4 = a+4.*h #we need to move four times the step size away from a.
    x5 = a+5.*h #we need to move five times the step size away from a.
    sum = 0.3*h*(11.*f(x1)-14.*f(x2)+26.*f(x3)-14.*f(x4)+11.*f(x5)) #formula
    return sum


#Problem 5.3.7

#let's test the open Newton-Cotes rules. 

def f(x): #to calculate the above function, we need an open interval (-1,1). There are asymptotes -1 and 1.
    return 1.0/((1.0-(x**2.0))**(0.5))

print nc_midpoint(f,-1,1,10) #returns 2.0
print nc_two_point(f,-1,1,10) #returns 2.12132010491
print nc_three_point(f,-1,1,10) #returns 2.41253476298
print nc_four_point(f,-1,1,10) #returns 2.46177011709
print nc_five_point(f,-1,1,10) #returns 2.58176125023

#The true result is 3.14159, or pi. You can see here that even with the five-point rule, we are still far away from
#the true result. This is due to the nature of the function, which as vertical asymptotes at x=-1 and x=1.
#Theoretically, adding additional points should improve the precision of the result, but this is impractical
#given other methods of integration, such as Gaussian quadrature or Richardson extrapolation. 



#Romberg Algorithm#

#The romberg algorithm produces a triangular array of numbers, all of which are numerical estimates
# of the definite integral of a function f(x) on the interval [a,b], where f(x) is continuous on [a,b].
#the romberg algorithm itself is an application of Richardson extrapolation to each of the iterative trapezoid approximations.
#This allows us to obtain a higher order extrapolation and thus a better result.


def romberg(f,a,b,n): #r is an array, f is a function on [a,b], n = number of rows.
    r = zeros((n,n), dtype='float') #declare an empty array of nxn dimension
    h = (b-a)
    r[0,0] = (h*0.5)*(f(a)+f(b)) #this works correctly.
    for i in range(1,n): #(1,n)?
        h = 0.5*h #takes previously calculated h and halves it.
        sum = 0
        for k in range(1,((2**i)),2):
            sum = sum + f(a+k*h)
        r[i,0] = 0.5*r[i-1,0]+(sum*h)
        for j in range(1,i+1):
            r[i,j] = (r[i,j-1]+((r[i,j-1]-r[i-1,j-1])/((4**j)-1)))
    return r
    
    
#Problem 5.2.2 and 5.2.3#

def f(x):
    return 4/(1+(x**2))

print romberg(f,0,1,6)

#let's test some other functions.

def g(x):
    return ((1+x)**(-1))

def h(x):
    return e**x

def j(x):
    return ((1+(x**2))**(-1))

print romberg(g,0,1,6) #correct result is: 0.69314
print romberg(h,0,1,6) #correct result is: 1.71828
print romberg(j,0,1,6) #correct result is: 0.78539

#let's try a bad function

def k(x):
    return sqrt(x)

print romberg(k,0,1,6) #correct value is 0.66667.

#k(x) is a bad function because the above algorithm can't handle estimation at the first endpoint, a=0.
#The result of that is that we need a far greater n to approach the solution than with other functions.
#increasing to n=20 allows us to approach the solution. This isn't ideal though... Increasing n even more
#slows down the algorithm immensely and, subsequently, a very large matrix is needed to store the result.


#trapezoid integration method
def trapint(f,a,b,n):
    stepsize = (b - a) / float(n)
    output = 0.5 *(f(a) + f(b))
    for i in range(1, n):+
        output += f(a + float(i) * stepsize)
        output *= stepsize
    return output


f = lambda x: math.e**(-1.*x**2.)

print trapint(f, 0., 1., 60)

f = lambda x: (1. + x**2.)**-1.

a = trapint(f, 0., 1., 2)

print a

print "actual value is atan(1), so the error is " + str(math.atan(1) - a)

print "according to the error formula, the maximum error for this fxn is " + str(1./24.)


#romberg algorithm
def romberg(f, a, b, n):
    h = b - a
    r = [[(h/2.)*(f(a) + f(b))]]
    for i in range(1, n+1):
        h /= 2.
        romsum = 0.
        for k in range(1, 2**i, 2):
            romsum += f(a + k * h)
        r.append([0.5*r[i-1][0] + romsum * h])
        for j in range(1,i + 1):
            r[i].append(r[i][j-1] + (r[i][j-1] - r[i-1][j-1])/(4.**j - 1.))
    return r

f = lambda x: math.e**(-1.*x**2.)

r = romberg(f, 0., 1., 5)

print r 


#Gaussian Quadrature
#Two point Gaussian Quadrature

def two_pt_gauss(f,a,b,n):
    h = (b-a)/n
    #print h
    sum = 0
    for i in range(n):
        x0 = a+(i*h) #starting at left end point, h represents step size.
        #print x0
        x1 = x0+(0.5*h)*(1-sqrt(1./3.))
        #print x1
        x2 = x0+(0.5*h)*(1+sqrt(1./3.))
        #print x2
        sum += ((f(x1)+f(x2))) #weights are 1.
        #print sum
    sum *= (0.5*h)
    return sum



#Three point Gaussian Quadrature#
def three_pt_gauss(f,a,b,n):
    h = (b-a)/n
    #print h
    sum = 0
    for i in range(n):
        x0 = a+(i*h) #starting at left end point, h represents step size.
        #print x0
        x1 = x0+(0.5*h)*(1-sqrt(3./5.))
        #print x1
        x2 = x0+(0.5*h)
        #print x2
        x3 = x0+(0.5*h)*(1+sqrt(3./5.))
        #print x3
        sum += (((5./9.)*f(x1))+((8./9.)*f(x2))+((5./9.)*f(x3))) #weights are 5/9, 8/9, and 5/9.
        #print sum
    sum *= (0.5*h)
    return sum


def f(x):
    return x**5

print two_pt_gauss(f,0.,1.,10) #result is 1.66665277778
print three_pt_gauss(f,0.,1.,10) #result is 1.66666666667

#You can see here that using Gaussian Quadrature feels like cheating, almost. With so little computation,
#we can approximate the area under the curve (i.e. integrate) with a, compared to other methods, very high
#amount of accuracy. Using n=10 and three_pt_gauss, we obtain the true result.


#Runge-Kutta#

#The Runge-Kutta method imitate the Taylor series method without requiring analytic differentiation of the
#original differential equation. Therefore, the below algorithms will accept any function, interval, and given point
#and will represent the solution to an ODE locally at the given point. This method works provided we now the
#value exactly at some arbitrary t. The downside of these methods is that they have to evaluate the function f several times,
#which can be very time consuming depending on the function.

#f denotes the function, x denotes the initial point, a and b denote the beginning and end of the interval, n denotes number of iterations.

#NOTE: if you wish to plot two functions, you must do so seperately.
#Plot one function, terminate the function, then plot the other.

# Runge-Kutta of order 2.
def runge_kutta_2(f, x, a, b, n):
    h = (b - a) / n  # step size
    t = a  # sets the initial t as the left-most point of the interval, a.
    for j in range(1, n + 1):
        k1 = h * f(t, x)
        k2 = h * f(t + h, x + k1)
        x += 0.5 * (k1 + k2)
        t = a + (j * h)
        print [j, t, x]


# Runge-Kutta method of order 4.
def runge_kutta_4(f, x, a, b, n):
    h = (b - a) / n  # step size
    t = a  # sets the initial t as the left-most point of the interval, a.
    for j in range(1, n + 1):
        k1 = h * f(t, x)
        k2 = h * f(t + 0.5 * h, x + 0.5 * k1)
        k3 = h * f(t + 0.5 * h, x + 0.5 * k2)
        k4 = h * f(t + h, x + k3)
        x += ((1. / 6.) * (k1 + 2. * k2 + 2. * k3 + k4))
        t = a + (j * h)
        print [j, t, x]


# Runge-Kutta method of order 4 that also plots.
def runge_kutta_4_plot(f, x, a, b, n):
    h = (b - a) / n  # step size
    t = a  # sets the initial t as the left-most point of the interval, a.
    jstore = []  # empty list to store j values. Not strictly necessary for plotting but may be useful.
    tstore = []  # empty list to store t values
    xstore = []  # empty list to store x values
    for j in range(1, n + 1):
        k1 = h * f(t, x)
        k2 = h * f(t + 0.5 * h, x + 0.5 * k1)
        k3 = h * f(t + 0.5 * h, x + 0.5 * k2)
        k4 = h * f(t + h, x + k3)
        x += ((1. / 6.) * (k1 + 2. * k2 + 2. * k3 + k4))
        t = a + (j * h)
        print [j, t, x]
        jstore.append(j)
        tstore.append(t)
        xstore.append(x)
    plt.plot(tstore, xstore)
    plt.xlabel('t')
    plt.ylabel('x')
    plt.show()

#For the Runge-Method of order 5 given below,
#the difference between the values of x(t+h) obtained from the 4th and 5th order procedures is an estimate of the
#local truncation error in the 4th order procedure. Therefore, evaluations give a 5th order approximation,
#together with an error estimate. This method is called the Runge-Kutta-Fehlberg method

#adaptive (enter 1 for yes, 0 for no) is a toggle to store error values needed in the adaptive Runge-Kutta method.

# Runge-Kutta method of order 5.
def runge_kutta_fehlberg(f, x, a, b, n, adaptive):
    h = (b - a) / n  # step size
    t = a  # sets the initial t as the left-most point of the interval, a.
    if adaptive == 1:
        erstore = [] #need a place to store error values for adaptive Runge-Kutta.
    #note that a2 == b2 == 0.
    #using long format to increase accuracy of estimation.
    c20 = long(0.25)
    c21 = long(0.25)
    c30 = long(0.375)
    c31 = long(0.09375)
    c32 = long(0.28125)
    c40 = long(12. / 13.)
    c41 = long(1932. / 2197.)
    c42 = long(7200. / 2197.)
    c43 = long(7296. / 2197.)
    c51 = long(439. / 216.)
    c52 = long(-8.)
    c53 = long(3680. / 513.)
    c54 = long(-845. / 4104.)
    c60 = long(0.5)
    c61 = long(-8. / 27.)
    c62 = long(2.)
    c63 = long(-3544. / 2565.)
    c64 = long(1859. / 4104.)
    c65 = long(-0.275)
    a1 = long(25. / 216.)
    a3 = long(1408. / 2565.)
    a4 = long(2197. / 4104.)
    a5 = long(-0.2)
    b1 = 16. / 135.
    b3 = 6656. / 12825.
    b4 = 28561. / 56430.
    b5 = -0.18
    b6 = 2. / 55.
    for j in range(1, n + 1):
        k1 = h * f(t, x)
        k2 = h * f(t + (c20 * h), x + (c21 * k1))
        k3 = h * f(t + (c30 * h), x + (c31 * k1) + (c32 * k2))
        k4 = h * f(t + (c40 * h), x + (c41 * k1) + (c42 * k2) + (c43 * k3))
        k5 = h * f(t + h, x + (c51 * k1) + (c52 * k2) + (c53 * k3) + (c54 * k4))
        k6 = h * f(t + (c60 * h), x + (c61 * k1) + (c62 * k2) + (c63 * k3) + (c64 * k4) + (c65 * k5))
        x4 = x + ((a1 * k1) + (a3 * k3) + (a4 * k4) + (a5 * k5))
        x += ((b1 * k1) + (b3 * k3) + (b4 * k4) + (b5 * k5) + (b6 * k6))
        t = a + (j * h)
        er = abs(x - x4)
        if adaptive == 1:
            erstore.append(er)
        else:
            print [j, t, x, er]
    if adaptive == 1:
        return erstore

#The Bogacki Shampine method is a Runge Kutta method of order three with four stages with the First
#Same As Last (FSAL) property, so that it uses approximately three function evaluations per step.
#It has an embedded second-order method which can be used to implement adaptive step size.

def runge_kutta_bogacki_shampine(f, x, a, b, n, adaptive):
    h = (b - a) / n  # step size
    t = a  # sets the initial t as the left-most point of the interval, a.
    if adaptive == 1:
        erstore = [] #need a place to store error values for adaptive Runge-Kutta.
    #note that b4 == c31 == 0
    #using long format to increase accuracy of estimation.
    c20 = long(0.5)
    c21 = long(0.5)
    c30 = long(0.75)
    c32 = long(0.75)
    c40 = long(1.)
    c41 = long(2. / 9.)
    c42 = long(1. / 3.)
    c43 = long(4. / 9.)
    a1 = long(2. / 9.)
    a2 = long(1. / 3.)
    a3 = long(4. / 9.)
    b1 = 7. / 24.
    b2 = 1. / 4.
    b3 = 1. / 3.
    b4 = 1. / 8.
    for j in range(1, n + 1):
        k1 = h * f(t, x)
        k2 = h * f(t + (c20 * h), x + (c21 * k1))
        k3 = h * f(t + (c30 * h), x + (c32 * k2))
        k4 = h * f(t + (c40 * h), x + (c41 * k1) + (c42 * k2) + (c43 * k3))
        x3 = x + ((a1 * k1) + (a2 * k2) + (a3 * k3))
        x += ((b1 * k1) + (b2 * k2) + (b3 * k3) + (b4 * k4))
        t = a + (j * h)
        er = abs(x - x3)
        if adaptive == 1:
            erstore.append(er)
        else:
            print [j, t, x, er]
    if adaptive == 1:
        return erstore

#The Cash-Karp Runge-Kutta method uses six function evaluations to calculate fourth and fifth order accurate solutions.
#The difference between these solutions is then taken to be the error of the fourth order solution. 

def runge_kutta_cash_karp(f, x, a, b, n, adaptive):
    h = (b - a) / n  # step size
    t = a  # sets the initial t as the left-most point of the interval, a.
    if adaptive == 1:
        erstore = [] #need a place to store error values for adaptive Runge-Kutta.
    #note that a2 == a5 == b2 == 0..
    #using long format to increase accuracy of estimation.
    c20 = long(1. / 5.)
    c21 = long(1. / 5.)
    c30 = long(3./ 10.)
    c31 = long(3./ 40.)
    c32 = long(9./ 40.)
    c40 = long(3. / 5.)
    c41 = long(3. / 10.)
    c42 = long(-9. / 10.)
    c43 = long(6. / 5.)
    c50 = long(1.)
    c51 = long(-11. / 54.)
    c52 = long(5. / 2.)
    c53 = long(-70. / 27.)
    c54 = long(35. / 27.)
    c60 = long(7. / 8.)
    c61 = long(1631. / 55296.)
    c62 = long(175. / 512.)
    c63 = long(575. / 13824.)
    c64 = long(44275. / 110592.)
    c65 = long(253. / 2096.)
    a1 = long(37. / 378.)
    a3 = long(250. / 621.)
    a4 = long(125. / 594.)
    a6 = long(512. / 1771.)
    b1 = 2825. / 27648.
    b3 = 18575. / 48384.
    b4 = 13525. / 55296.
    b5 = 277. / 14336.
    b6 = 1. / 4.
    for j in range(1, n + 1):
        k1 = h * f(t, x)
        k2 = h * f(t + (c20 * h), x + (c21 * k1))
        k3 = h * f(t + (c30 * h), x + (c31 * k1) + (c32 * k2))
        k4 = h * f(t + (c40 * h), x + (c41 * k1) + (c42 * k2) + (c43 * k3))
        k5 = h * f(t + (c50 * h), x + (c51 * k1) + (c52 * k2) + (c53 * k3) + (c54 * k4))
        k6 = h * f(t + (c60 * h), x + (c61 * k1) + (c62 * k2) + (c63 * k3) + (c64 * k4) + (c65 * k5))
        x4 = x + ((a1 * k1) + (a3 * k3) + (a4 * k4) + (a6 * k6))
        x += (b1 * k1) + (b3 * k3) + (b4 * k4) + (b5 * k5) + (b6 * k6)
        t = a + (j * h)
        er = abs(x - x4)
        if adaptive == 1:
            erstore.append(er)
        else:
            print [j, t, x, er]
    if adaptive == 1:
        return erstore

#The Dormand Prince method is another type of ODE solver in the Runge-Kutta family.
#It uses six evaluations to calculate fourth and fifth order accurate solutions.
#The Dormand Pricnce has seven stages, but it uses only six function evaluations per step because
#it has the FSAL (First Same As Last) property: the last stage is evaluated at the same point as the first stage
#of the next step. Dormand and Prince chose the coefficients of their method to minimize the error of the fifth-order
#solution. This is the main difference with the Fehlberg method, which was constructed so that the fourth-order
#solution has small error. For this reason, the Dormand Prince method is more suitable when the higher-order
#solution is used to continue the integration.

def runge_kutta_dormand_prince(f, x, a, b, n, adaptive):
    h = (b - a) / n  # step size
    t = a  # sets the initial t as the left-most point of the interval, a.
    if adaptive == 1:
        erstore = [] #need a place to store error values for adaptive Runge-Kutta.
    #note that a2 == a6 == b2 == b6 == c72 == 0.
    #using long format to increase accuracy of estimation.
    c20 = long(1. / 5.)
    c21 = long(1. / 5.)
    c30 = long(3. / 10.)
    c31 = long(3. / 40.)
    c32 = long(9. / 40.)
    c40 = long(4. / 5.)
    c41 = long(44. / 45.)
    c42 = long(-56. / 15.)
    c43 = long(32. / 9.)
    c50 = long(8. / 9.)
    c51 = long(19372. / 6561.)
    c52 = long(-25360. / 2187.)
    c53 = long(64448. / 6561.)
    c54 = long(-212. / 729.)
    c60 = long(1.)
    c61 = long(9017. / 3168.)
    c62 = long(-355. / 33.)
    c63 = long(46732. / 5247.)
    c64 = long(49. / 176.)
    c65 = long(-5103. / 18656.)
    c70 = long(1.)
    c71 = long(35. / 384.)
    c73 = long(500. / 1113.)
    c74 = long(125. / 132.)
    c75 = long(-2187. / 6784.)
    c76 = long(11. / 84.)
    a1 = long(35. / 384.)
    a3 = long(500. / 1113.)
    a4 = long(125. / 192.)
    a5 = long(-2187. / 6784.)
    b1 = 5179. / 57600.
    b3 = 7571. / 16695.
    b4 = 393. / 640.
    b5 = -92097. / 339200.
    b6 = 187. / 2100.
    b7 = 1. / 40.
    for j in range(1, n + 1):
        k1 = h * f(t, x)
        k2 = h * f(t + (c20 * h), x + (c21 * k1))
        k3 = h * f(t + (c30 * h), x + (c31 * k1) + (c32 * k2))
        k4 = h * f(t + (c40 * h), x + (c41 * k1) + (c42 * k2) + (c43 * k3))
        k5 = h * f(t + (c50 * h), x + (c51 * k1) + (c52 * k2) + (c53 * k3) + (c54 * k4))
        k6 = h * f(t + (c60 * h), x + (c61 * k1) + (c62 * k2) + (c63 * k3) + (c64 * k4) + (c65 * k5))
        k7 = h * f(t + (c70 * h), x + (c71 * k1) + (c73 * k3) + (c74 * k4) + (c75 * k5) + (c76 * k6))
        x5 = x + ((a1 * k1) + (a3 * k3) + (a4 * k4) + (a5 * k5))
        x += ((b1 * k1) + (b3 * k3) + (b4 * k4) + (b5 * k5) + (b6 * k6) + (b7 * k7))
        t = a + (j * h)
        er = abs(x - x5)
        if adaptive == 1:
            erstore.append(er)
        else:
            print [j, t, x, er]
    if adaptive == 1:
        return erstore

#The error estimate 'er' calculated and stored in some of the above methods
#can tell us when to adjust the step size to control for the single-step error.
#This fact, together with the 5th order approximation, results in an adaptive procedure
#that is very accurate. emin and emax are the lower and upper bounds on the allowable error estimate.
#hmin and hmax are are bounds on the step size h.
#eflag is an error flag that returns 0 or 1 depending on if there is a successful march from a to b or if max. n is reached.
#method is a toggle to select between different Runge-Kutta methods, as follows.

#1: Runge-Kutta-Fehlberg
#2: Runge-Kutta Bogacki-Shampine
#3: Runge-Kutta Cash-Karp
#4: Runge-Kutta Dormand-Prince

def runge_kutta_adaptive(f, x, a, b, n, emin, emax, hmin, hmax, method):
    erstore = []
    if method == 1:
        erstore = array(runge_kutta_fehlberg(f, x, a, b, n, 1)) #store error values for use in adaptive method.
    if method == 2:
        erstore = array(runge_kutta_bogacki_shampine(f, x, a, b, n, 1)) #store error values for use in adaptive method.
    if method == 3:
        erstore = array(runge_kutta_cash_karp(f, x, a, b, n, 1)) #store error values for use in adaptive method.
    if method == 4:
        erstore = array(runge_kutta_dormand_prince(f, x, a, b, n, 1)) #store error values for use in adaptive method.
    h = (b - a) / n
    t = a
    eflag = 1
    k = 0
    sig = 0.5 * (10.0 ** -5.)
    while k <= n:
        #maybe add for loop here for erstore at the bottom of the function.
        k = k + 1
        if abs(h) < hmin:
            h = copysign(hmin,(h)) #return h with sign(h)*hmin
        if abs(h) > hmax:
            h = copysign(hmax,(h)) #return h with sign(h)*hmax
        d = abs(b - a)
        if d <= abs(h):
            eflag = 0
            if d <= sig * max(abs(b), abs(a)):
                break
            h = copysign(d,(h)) #return h with sign(h)*d
        xsave = x
        tsave = t
        if method == 1:
            print runge_kutta_fehlberg(f, x, a, b, n, 0) #note: no need to return error values for printing results
        if method == 2:
            print runge_kutta_bogacki_shampine(f, x, a, b, n, 0)
        if method == 3:
            print runge_kutta_cash_karp(f, x, a, b, n, 0)
        if method == 4:
            print runge_kutta_dormand_prince(f, x, a, b, n, 0)
        if eflag == 0:
            break
        for i in range(1, n + 1):
            if i in erstore < emin:
                h = 2. * h
            if i in erstore > emax:
                h = 0.5 * h
                x = xsave
                t = tsave
                k = k - 1
                
# test function given in book.
def f(t, x):
    return 2. + ((x - t - 1.) ** 2.)

#another test function
def g(t,x):
    return 3. + (5 * sin(t)) + (0.2 * x)

#tests using f(t,x):
print runge_kutta_2(f, 2., 1., 1.5625, 72)
print runge_kutta_4(f, 2., 1., 1.5625, 72) #returns 3.192937673837072
print runge_kutta_4_plot(f, 2., 1., 1.5625, 72) #same as above
print runge_kutta_fehlberg(f, 2., 1., 1.5625, 72, 0) #returns 3.21488255297631
print runge_kutta_bogacki_shampine(f, 2., 1., 1.5625, 72, 0) #returns 3.1907582960278154
print runge_kutta_cash_karp(f, 2., 1., 1.5625, 72, 0 #returns 3.192865860211696)
print runge_kutta_dormand_prince(f, 2., 1., 1.5625, 72, 0) #returns 3.19079031694774
print runge_kutta_adaptive(f, 2., 1., 1.5625, 72, 10 ** -8, 10 ** -5, 10 ** -6, 1.0, 4) #returns 3.19079031694774

#tests using g(t,x):
print runge_kutta_2(g, 0., 0., 10., 1000)
print runge_kutta_4(g, 0., 0., 10., 1000) #returns 135.91724460992168.
print runge_kutta_4_plot(g, 0., 0., 10., 1000) #same as above
print runge_kutta_fehlberg(g, 0., 0., 10., 100000, 0) #returns 135.93269120019045.
print runge_kutta_bogacki_shampine(g, 0., 0., 10., 100000, 0) #returns 135.91374364421674
print runge_kutta_cash_karp(g, 0., 0., 10., 100000, 0) #returns 135.91543072388197
print runge_kutta_dormand_prince(g, 0., 0., 10., 100000, 0) #returns 135.91373784394202
print runge_kutta_adaptive(g, 0., 0., 10., 1000, 10 ** -8, 10 ** -5, 10 ** -6, 1.0, 4) #returns 135.56713499140957

def g(t, x):
    return -2. * t

#print runge_kutta_4(g,1.,0.,e,100) #we obtain -6.3890560989306495, correct result is = -6.3890560989
print runge_kutta_4_plot(g, 1., 0., e, 100)

def f(t, x):  # the function sin(t)/t is not defined at t==0, so we must specify a value at t==0.
    if t == 0:
        return 1.
    else:
        return sin(t) / t

print runge_kutta_4(f, 0., 0., 1., 100)  #(1-0)/100 = 0.01 so step size is as the problem specifies. Result = 0.9460830703677978. Correct result is 0.9460830703.

def f(t, x):
    return (2. / sqrt(pi)) * (e ** (-t ** 2.))

# let's approximate this by g(t0, t2, n) and see how close it is to our RK4 approximation.

def g(t0, t2, n):
    h = (t2 - t0) / n
    for j in range(1, n + 1):
        t = t0 + (j * h)
        a = 0.3084284
        b = -0.0849713
        c = 0.6627698
        y = (1. + (0.47047 * t)) ** (-1.)
        x = 1. - ((a * y) + (b * (y ** 2.)) + (c * (y ** 3.))) * ((2. / sqrt(pi)) * (e ** (-t ** 2.)))
        print [j, t, x]


print runge_kutta_4(f, 0., 0., 2., 100)  #result is 0.9953222649730256
print g(0., 2., 100)  #result is 0.9953087431727834.

#so our RK4(f(t,x)) and g(t0,t2,n) evaluations differ by 0.0001352 at 2.
#The g(t0,t2,n) function is a reasonably good approximation of our complicated f(t,x) function.

#A slightly different Runge-Kutta 4 algorithm. 
#here's my rk4 function, it seems to work well
def rk4(yp, y, a, b, n):
    h = (b - a) / float(n)
    xlist = [a]
    ylist = [y]
    for i in range(1, n + 1):
        k1 = h * yp(xlist[i-1], ylist[i-1])
        k2 = h * yp(xlist[i-1] + .5*h, ylist[i-1] + .5*k1)
        k3 = h * yp(xlist[i-1] + .5*h, ylist[i-1] + .5*k2)
        k4 = h * yp(xlist[i-1] + h, ylist[i-1] + k3)
        ylist.append(ylist[i-1] + (1. / 6.)*(k1 + 2. * (k2 + k3) + k4))
        xlist.append(a + i * h)
    return xlist, ylist
    
yp = lambda x, y: 2. + (y - x - 1) ** 2.
y = 2.
a = 1.
b = 1.5625
n = 72

x, y = rk4(yp, y, a, b, n)

plt.plot(x, y)




#convert polar coordinates to Cartesian coordinates#

# Ask the user for the values of r and theta
r = float(input("Enter r: "))
d = float(input("Enter theta in degrees: "))

# Convert the angle to radians
theta = d*pi/180

# Calculate the equivalent Cartesian coordinates
x = r*cos(theta)
y = r*sin(theta)

# Print out the results
print("x = ",x,", y = ",y)


#Print a Fibonacci sequence#
f1,f2 = 1,1
while f2<1000:
    print(f2)
    f1,f2 = f2,f1+f2


#Jacobian Matrix Evaluation#

#given two (or three) functions with respect to x and y (and z), calculates the determinant of the Jacobian matrix at a given point (xi,yi).
#useful in characterizing local behavior of nonlinear system about an equilibrium point, i.e. linearization.

#fx,fy,gx,gy are functions f and g with respect to x and y
#xi,yi are points of which to evaluate the Jacobian

def jacobian_eval_R2(fx,fy,gx,gy,xi,yi):
    a = fx(xi,yi)
    b = fy(xi,yi)
    c = gx(xi,yi)
    d = gy(xi,yi)
    Jcomplete = [[a,b],[c,d]]
    return linalg.det(Jcomplete)

def jacobian_eval_R3(fx,fy,fz,gx,gy,gz,hx,hy,hz,xi,yi,zi):
    a = fx(xi,yi,zi)
    b = fy(xi,yi,zi)
    c = fz(xi,yi,zi)
    d = gx(xi,yi,zi)
    e = gy(xi,yi,zi)
    f = gz(xi,yi,zi)
    g = hx(xi,yi,zi)
    h = hy(xi,yi,zi)
    i = hz(xi,yi,zi)
    Jcomplete = [[a,b,c],[d,e,f],[g,h,i]]
    return linalg.det(Jcomplete)

#test functions for R2
def f(x,y):
    return 2*x*(y**3)

def g(x,y):
    return y*(x**2)

def fx(x,y):
    return 2*(y**3)

def fy(x,y):
    return 6*x*(y**2)

def gx(x,y):
    return 2*x*y

def gy(x,y):
    return x**2

print jacobian_eval_R2(fx,fy,gx,gy,1.,1.)


#test functions for R3
def f3(x,y,z):
    return 2*x*(y**3)*z

def g3(x,y,z):
    return y*(x**2)*(z**4)

def h3(x,y,z):
    return (y**4)*(z**3)

def fx3(x,y,z):
    return 2*(y**3)*z

def fy3(x,y,z):
    return 6*x*(y**2)*z

def fz3(x,y,z):
    return 2*x*(y**3)

def gx3(x,y,z):
    return 2*x*y*(z**4)

def gy3(x,y,z):
    return x**2*(z**4)

def gz3(x,y,z):
    return 4*(z**3)*y*(x**2)

def hx3(x,y,z):
    return 0

def hy3(x,y,z):
    return 4*(y**3)*(z**3)

def hz3(x,y,z):
    return 3*(z**2)*(y**4)

print jacobian_eval_R3(fx3,fy3,fz3,hx3,gy3,gz3,hx3,hy3,hz3,1.,1.,1.)

#The Wronskian at a given point#
#Ideally, we would evaluate at a range of points. Future version may incorporate this.

#need to define f, f', g, g', and point (x,y).
def wronskian_eval_R2(f,fprime,g,gprime,x,y): #wronskian for a system of two equations.
    a = f(x,y)
    b = fprime(x,y)
    c = g(x,y)
    d = gprime(x,y)
    w = [[a,b],[c,d]]
    if linalg.det(w) == 0:
        print "linearly dependent"
    else:
        print "linearly independent"
    return linalg.det(w)


def f(x,y):
    return x**2

def fprime(x,y):
    return 2*x*(y**3)

def g(x,y):
    return x**3

def gprime(x,y):
    return 3*(x**2)*y

print wronskian_eval_R2(f,fprime,g,gprime, 1., 1.) #this is the Wronskian only at one point, so apply Abel's theorem.
