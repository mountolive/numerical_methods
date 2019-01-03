from __future__ import division
from math import sqrt


class FirstAttemptMuller:
    ''' Adapted from
    https://github.com/apauley/numerical-analysis/blob/master/Chapter1/Python/muller.py'''

    @staticmethod
    def swap_points(x):
        s = x
        s.sort()
        f = s[1]
        sn = s[2]
        t = s[0]
        s[0] = f
        s[1] = sn
        s[2] = t
        return s

    def mullers_method(self, func, xl, xm, xr, tolerance=1e-2, n=100):
        x = [xl, xm, xr]
        i = 0
        while abs(x[0] - x[2]) > tolerance:
            i += 1
            x = self.swap_points(x)
            y = func(x[0]), func(x[1]), func(x[2])
            h1 = x[1]-x[0]
            h2 = x[0]-x[2]
            lam = h2/h1
            c = y[0]
            a = (lam*y[1] - y[0]*(1.0 + lam)+y[2])/(lam*h1**2.0*(1+lam))
            b = (y[1] - y[0] - a*(h1**2.0))/h1
            if b > 0:
                root = x[0] - ((2.0*c)/(b + (b**2 - 4.0*a*c)**0.5))
            else:
                root = x[0] - ((2.0*c)/(b - (b**2 - 4.0*a*c)**0.5))
            if abs(func(root)) > x[0]:
                x = [x[1], x[0], root]
            else:
                x = [x[2], x[0], root]
            x = self.swap_points(x)
            if i > n:
                print("Max iterations permitted surpassed, aborting Muller's")
                break
        print("Iterations for convergence (if > n, the result is not the root) #: %d" % i)
        return x[0]


class SecondMuller:
    '''From: http://www2.gsu.edu/~matrhc/muller.py'''

    def muller(self, f, p0, p1, p2, tol=1e-2, n=100):
        h1 = p1 - p0
        h2 = p2 - p1
        f_p1 = f(p1)
        f_p2 = f(p2)
        d1 = (f_p1 - f(p0)) / h1
        d2 = (f_p2 - f_p1) / h2
        d = (d2 - d1) / (h2 + h1)
        i = 1
        p = 0.0
        while i <= n:
            b = d2 + h2 * d
            # sqrt function knows to use complex numbers
            # if we make sure the argument is complex
            # (by adding 0j in case argument starts as real)
            # D is the discriminant
            D = sqrt(b * b - 4 * f_p2 * d + 0j)
            if abs(b - D) < abs(b + D):
                # ensure we don't subtract similar values
                E = b + D
            else:
                E = b - D
            h = -2 * f_p2 / E
            p = p2 + h
            if abs(h) < tol:
                print("Convergence reached at %d" % i)
                if not p.imag:
                    p = p.real
                return p
            p0 = p1
            p1 = p2
            p2 = p
            h1 = p1 - p0
            h2 = p2 - p1
            f_p1 = f(p1)
            f_p2 = f(p2)
            d1 = (f_p1 - f(p0)) / h1
            d2 = (f_p2 - f_p1) / h2
            d = (d2 - d1) / (h2 + h1)
            i += 1
        print("Reached maximum number of iterations")
        return p
