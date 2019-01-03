import csv
import math


class NumericalMethod:
    '''
    This module defines the bisection and Muller's methods for finding
    roots of a given function. It also prints the results into a csv
    '''

    def __init__(self):
        self.iterations = []

    def bisection(self, func, xl=0.1, xr=56.0, tolerance=1e-2, n=1000):
        self._clean_list()
        err = 1
        i = 0
        while err > tolerance:
            func_left, func_right = func(xl), func(xr)
            row = [xl, xr, func_left, func_right]
            mid = (xr + xl) / 2.0
            func_mid = func(mid)
            row += [mid, func_mid]
            if self.same_sign(func_left, func_mid):
                xl = mid
            else:
                xr = mid
            i += 1
            err = abs(xl - xr)
            row.append(err)
            self.iterations.append(row)
            i += 1
            if i > n:
                print('Max number of iteration hit. Method unsuccessful')
                break
        return xl

    def muller(self, func, xl=0.1, xr=56.0, tolerance=1e-2, n=1000):
        # TODO: needs to be fixed
        self._clean_list()
        err = 1
        i = 0
        while err > tolerance:
            xm = (xl + xr) / 2.0
            a, b, c = self._critical_points_muller(func, xl, xm, xr)
            row = [xr, xm, xl, a, b, c, err]
            x1, x2 = self._quadratic_roots(a, b, c)
            complex_1 = type(x1) is complex
            complex_2 = type(x2) is complex
            if complex_1 and complex_2:
                self.iterations.append(row)
                print('Stopped by complex root!')
                break
            elif not complex_1 and abs(func(x1)) < tolerance:
                self._refactor_list(x1, row, a, b, c, tolerance)
                break
            elif not complex_2 and abs(func(x2)) < tolerance:
                self._refactor_list(x1, row, a, b, c, tolerance)
                break
            elif complex_1:
                root = x2
            elif complex_2:
                root = x1
            else:
                root = min(x1, x2)
            xl, xm, xr = self._new_roots(root, xl, xm, xr)
            print(xl, xm, xr)
            err = abs(xl - xr)
            row.append(err)
            self.iterations.append(row)
            i += 1
            if i > n:
                print('Max number of iteration hit. Method unsuccessful')
                break
        return xl

    def regula_falsi(self, func, xl=0.1, xr=56.0, tolerance=1e-2, n=1000):
        err = 1
        i = 0
        while err > tolerance:
            yr = func(xr)
            yl = func(xl)
            row = [xl, xr, yl, yr]
            xprime = self._regula_root(func, xl, xr)
            yprime = func(xprime)
            if not self.same_sign(yprime, yr):
                xl = xprime
            elif not self.same_sign(yprime, yl):
                xr = xprime
            else:
                print('Found same sign interval, breaking')
                break
            err = abs(xr - xl)
            row.append(err)
            self.iterations.append(row)
            i += 1
            if i > n:
                print('Max number of iteration hit. Method unsuccessful')
                break
        return xl

    def print_csv(self, name, method='bisection'):
        with open(name, 'w') as csvfile:
            writer = csv.writer(csvfile)
            if method == 'bisection':
                writer.writerow(['iteration', 'Xleft', 'Xright', 'f(Xl)', 'f(Xr)', 'f(mid)', 'Xmid'
                                 'tolerance'])
            elif method == 'muller':
                writer.writerow(['iteration', 'Xleft', 'Xmid', 'Xright', 'a', 'b', 'c', 'tolerance'])
            elif method == 'regula':
                writer.writerow(['iteration', 'Xleft', 'Xright', 'f(Xl)', 'f(Xr)', 'tolerance'])
            else:
                raise Exception('Method passed not available')
            for i, row in enumerate(self.iterations):
                row = [i] + row
                writer.writerow(row)

    def formated_print(self, method='bisection'):
        header = ' {:^19s} |'
        rows = ' {:^19} |'
        if method == 'bisection':
            header *= 8
            rows *= 8
            print(header.format('iteration', 'Xleft', 'Xright', 'f(Xl)', 'f(Xr)', 'f(mid)', 'Xmid',
                                'tolerance'))
        elif method == 'regula':
            header *= 6
            rows *= 6
            print(header.format('iteration', 'Xleft', 'Xright', 'f(Xl)', 'f(Xr)', 'tolerance'))
        else:
            header *= 8
            rows *= 8
            print(header.format('iteration', 'Xleft', 'Xmid', 'Xright', 'a', 'b', 'c', 'tolerance'))
        for i, row in enumerate(self.iterations):
            print(rows.format(*([i] + row)))

    @staticmethod
    def same_sign(a, b):
        return a * b > 0

    @staticmethod
    def _regula_root(func, xl, xr):
        yr = func(xr)
        yl = func(xl)
        return (xl * yr - xr * yl) / (yr - yl)

    @staticmethod
    def _critical_points_muller(func, xl, xm, xr):
        ul = xl - xm
        ur = xr - xm
        fr = func(xr)
        fm = func(xm)
        fl = func(xl)
        a = ((fr - fm) * ul - (fl - fm) * ur) / (ur**2 * ul - ur * ul**2)
        b = ((fl - fm) * ur**2 - (fr - fm) * ul) / (ur**2 * ul - ur * ul**2)
        c = fm
        return a, b, c

    @staticmethod
    def _quadratic_roots(a, b, c):
        root_1 = (-b + (b**2 - 4*a*c)**0.5) / (2 * a)
        root_2 = (-b - (b**2 - 4*a*c)**0.5) / (2 * a)
        return root_1, root_2

    @staticmethod
    def _new_roots(root, xl, xm, xr):
        roots = [xl, xm, xr]
        roots.sort(key=lambda x: abs(root - x))
        roots = (roots[:-1] + [root])
        roots.sort()
        return (elem for elem in roots)

    def _refactor_list(self, root, row, a, b, c, tolerance):
        self.iterations.append(row)
        xl = xr = xm = root
        row = [xr, xm, xl, a, b, c, tolerance]
        self.iterations.append(row)

    def _clean_list(self):
        self.iterations = []


def minimum_angle_example(x):
    y = 180.0 - 123.0 - x
    return 3 * (math.cos(x) / math.sin(x)**2) - 2 * ((math.cos(y)) / math.sin(y)**2)


if __name__ == '__main__':
    solver = NumericalMethod()
    solver.bisection(minimum_angle_example)
    solver.formated_print()
    solver.bisection(minimum_angle_example, tolerance=1e-6)
    solver.formated_print()
    solver.muller(minimum_angle_example)
    solver.print_csv('muller.csv', method='muller')
    print("Done!")
