import csv
import math


class NumericalMethod:
    '''
    This module defines the bisection and Muller's methods for finding
    roots of a given function. It also prints the results into a csv
    '''

    def __init__(self):
        self.csv_rows = []

    def bisection(self, func, xl=0.1, xr=56.0, tolerance=1e-2):
        self._clean_list()
        err = 1
        i = 0
        while err > tolerance:
            func_left, func_right = func(xl), func(xr)
            row = [xl, xr, func_left, func_right]
            mid = (xr + xl) / 2.0
            func_mid = func(mid)
            row += [mid, func_mid]
            if self._same_sign(func_left, func_mid):
                xl = mid
            else:
                xr = mid
            i += 1
            err = abs(xl - xr)
            row.append(err)
            self.csv_rows.append(row)

    def muller(self, func, xl=0.1, xr=56.0, tolerance=1e-2):
        self._clean_list()
        err = 1
        while err > tolerance:
            xm = (xl + xr) / 2.0
            a, b, c = self._critical_points_muller(func, xl, xm, xr)
            row = [xr, xm, xl, a, b, c, err]
            x1, x2 = self._quadratic_roots(a, b, c)
            complex_1 = type(x1) is complex
            complex_2 = type(x2) is complex
            if complex_1 and complex_2:
                self.csv_rows.append(row)
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
            xl, xr = self._new_roots(root, xl, xm, xr)
            err = abs(xl - xr)
            row.append(err)
            self.csv_rows.append(row)

    def print_csv(self, name, method='bisection'):
        with open(name, 'w') as csvfile:
            writer = csv.writer(csvfile)
            if method == 'bisection':
                writer.writerow(['iteration', 'Xleft', 'Xright', 'f(Xl)', 'f(Xr)', 'f(mid)',
                                 'tolerance'])
            elif method == 'muller':
                writer.writerow(['iteration', 'Xleft', 'Xmid', 'Xright', 'a', 'b', 'c', 'tolerance'])
            else:
                raise Exception('Method passed not available')
            for i, row in enumerate(self.csv_rows):
                row = [i] + row
                writer.writerow(row)

    def _critical_points_muller(self, func, xl, xm, xr):
        ul = xl - xm
        ur = xr - xm
        fr = func(xr)
        fm = func(xm)
        fl = func(xl)
        a = ((fr - fm) * ul - (fl - fm) * ur) / (ur**2 * ul - ur * ul**2)
        b = ((fl - fm) * ur**2 - (fr - fm) * ul) / (ur**2 * ul - ur * ul**2)
        c = fm
        return a, b, c

    def _quadratic_roots(self, a, b, c):
        root_1 = (-b + (b * b - 4*a*c)**(1/2)) / (2 * a)
        root_2 = (-b - (b * b - 4*a*c)**(1/2)) / (2 * a)
        return root_1, root_2

    def _new_roots(self, root, xl, xm, xr):
        order = sorted([abs(root - xl), abs(root - xm), abs(root - xr)])[:-1]
        return (elem for elem in order)

    def _refactor_list(self, root, row, a, b, c, tolerance):
        self.csv_rows.append(row)
        xl = xr = xm = root
        row = [xr, xm, xl, a, b, c, tolerance]
        self.csv_rows.append(row)

    def _clean_list(self):
        self.csv_rows = []

    def _same_sign(self, a, b):
        return a * b > 0


def example_function(x):
    y = 180.0 - 123.0 - x
    return 3 * (math.cos(x) / math.sin(x)**2) - 2 * ((math.cos(y)) / math.sin(y)**2)


if __name__ == '__main__':
    solver = NumericalMethod()
    solver.bisection(example_function)
    solver.print_csv('biseccion.csv')
    solver.bisection(example_function, tolerance=1e-6)
    solver.print_csv('biseccion-6dec.csv')
    solver.muller(example_function)
    solver.print_csv('muller.csv', method='muller')
    print("Done!")
