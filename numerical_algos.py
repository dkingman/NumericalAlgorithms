# -*- coding: utf-8 -*-
"""
Created on Tue Jan 29 23:14:17 2013

@author: David Kingman
"""
import math
import numpy as np


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


def newton(f, f_prime, x, delta=1e-3, ap=1e-6, max_steps=25):
    fx = f(x)
    for n in range(1, max_steps):
        fp = f_prime(x)
        if abs(fp) < delta:
            print 'Small derivative'
            return [n, x, fx]
        d = fx / fp
        x = x - d
        #print [n,x,fx]
        fx = f(x)
        if abs(d) < ap:
            print 'Convergence'
            return [n, x, fx]


def func_3(x):
    return 2.0 * x * (1 - x ** 2 + x) * math.log(x) - x * x + 1


def func_3_prime(x):
    return -2.0 * (-1 + x) * (1 + x + (1 + 3 * x) * math.log(x))


print newton(func_3, func_3_prime, 0.9)


def newton_8(f, f_prime, x, delta=1e-6, ap=1e-6, max_steps=8):
    fx = f(x)
    for n in range(1, max_steps):
        fp = f_prime(x)
        if abs(fp) < delta:
            print 'Small derivative'
            return [n, x, fx]
        d = fx / fp
        x = x - d
        print [n, x, fx]
        fx = f(x)
        if abs(d) < ap:
            print 'Convergence'
            return [n, x, fx]


def func_4(x):
    return x ** 3 - math.sin(x) + 7


def func_4_prime(x):
    return x ** 2 - math.cos(x)


def func_5(x):
    return math.sin(x) - 1 + x


def func_5_prime(x):
    return math.cos(x) + 1


newton_8(func_4, func_4_prime, .5)
newton_8(func_5, func_5_prime, .5)


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


def func_6(x):
    return x * (x * x - 3.0) + 1


def func_6_prime(x):
    return 3.0 * x * x - 3


def func_7(x):
    return x ** 3 - 2 * math.cos(math.pi / 2 - x)


def func_7_prime(x):
    return -1.0 + x * x / 2 * (7 - x * x / 12)


def newton_2(f, f_prime, x, delta=1e-3, ap=1e-6, max_steps=25):
    fx = f(x)
    for n in range(1, max_steps):
        print 'n=', n, ' x=', x, ' fx=', fx
        fp = f_prime(x)
        if abs(fp) < delta:
            print 'Small derivative'
            return [n, x, fx]
        d = fx / fp
        x = x - d
        fx = f(x)
        if abs(d) < ap:
            print 'Convergence'
            return [n, x, fx]


a = newton_2(func_6, func_6_prime, 2.0)[1]
secant(func_6, a, 2.0)
a = newton_2(func_7, func_7_prime, 0.5)[1]
secant(func_7, a, 2.)


def secant_2(a, b, ap=1e-6, max_steps=25):
    for k in range(1, 11):
        f = lambda x: 2 * math.e ** (-k) * x + 1 - 3 * math.e ** (-k * x)
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
            #print 'n=',n,' a=',a,' fa=',fa
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
                print 'k=', k, ' n=', n, ' a=', a, ' fa=', fa
                break
            a = a - d
            fa = f(a)
            print 'k=', k, ' n=', n, ' a=', a, ' fa=', fa
            break


secant_2(0, 1)


def newton(f, f_prime, x, delta=1e-6, ap=1e-6, max_steps=25):
    fx = f(x)
    for n in range(1, max_steps):
        fp = f_prime(x)
        #print fx, fp
        if abs(fp) < delta:
            print 'Small derivative'
            return [n, x, fx]
        d = fx / fp
        x = x - d
        #print [n,x,fx]
        fx = f(x)
        if abs(d) < ap:
            print 'Convergence'
            return [n, x, fx]


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