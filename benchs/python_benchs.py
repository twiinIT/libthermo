import timeit
import numpy as np
from pythermo import IdealGas, RealGas
import json

ideal_g = IdealGas(287, 1004)
real_g = RealGas(287.05)

def bench_numpy(g, size, repeat, number):
    N = int(size)
    t1 = np.repeat(273., N)
    t2 = np.repeat(500., N)
    p1 = np.repeat(101325., N)
    p2 = np.repeat(5e5, N)
    eff = np.repeat(0.8, N)
    res = np.ones_like(t1)

    results = dict()

    print(f"Benchmark pythermo - RealGas - {N} elements")

    results["size"] = [N,]

    def bench(name, f):
        # run the bench and get exec time in nanosecond / element
        r = np.array(timeit.repeat(f, repeat=repeat, number=number)) * (1e9 / N / number)
        print(name, np.mean(r), np.std(r))
        results[name] = [np.mean(r),]
    
    bench("cp", lambda: g.cp(t1))
    bench("gamma", lambda: g.gamma(t1))
    bench("h", lambda: g.h(t1))
    bench("phi", lambda: g.phi(t1))
    bench("pr", lambda: g.pr(t1, t2, eff))
    bench("eff_poly", lambda: g.eff_poly(p1, t1, p2, t2))

    return results


def warm_up():
    t = np.repeat(273., 1e6)
    timeit.repeat(lambda: real_g.cp(t), repeat=10, number=500)

warm_up()

def merge_dicts(a, b, path=None):
    if path is None:
        path = []
    for key in b:
        if key in a:
            if isinstance(a[key], dict) and isinstance(b[key], dict):
                merge_dicts(a[key], b[key], path + [str(key)])
                continue

            if isinstance(a[key], list) and isinstance(b[key], list):
                a[key] = a[key] + b[key]
            elif a[key] == b[key]:
                pass  # same leaf value
            else:
                a[key] = b[key]
        else:
            a[key] = b[key]
    return a


real_results = dict()
# first are slow, last are large
#for s, n in ((1e1, 1000), (1e2, 1000), (1e3, 100), (1e4, 10000), (1e5, 1000), (1e6, 1000), (1e7, 100), (1e8, 10)):
for s, n in ((1e1, 1), (1e2, 1), (1e3, 1), (1e4, 1), (1e5, 100), (2.5e5, 100), (5e5, 100), (8e5, 100), (1e6, 100), (2e6, 1), (5e6, 1), (8e6, 1), (1e7, 1), (1e8, 1), (2e8, 1)):
    # merge_dicts(real_results, bench_numpy(real_g, s, repeat=10, number=n))
    pass

with open('benchs/benchs_pythermo_results.json', 'w') as f:
    json.dump(real_results, f, indent=2)


class NumpyRealGas:

    cp_coeff = [ 2.35822e-20,  -1.79613e-16, 4.70972e-13, -3.3094e-10,
                -6.27984e-07, 0.00123785,   -0.521742,   1068.43 ]
    h_coeff = [ 2.35822e-20 / 8., -1.79613e-16 / 7., 4.70972e-13 / 6.,
                -3.3094e-10 / 5., -6.27984e-07 / 4., 0.00123785 / 3.,
                -0.521742 / 2.,   1068.43,           0. ]
    phi_coeff = [ 2.35822e-20 / 7., -1.79613e-16 / 6.,
                4.70972e-13 / 5.,
                -3.3094e-10 / 4.,
                -6.27984e-07 / 3.,
                0.00123785 / 2.,
                -0.521742,
                0.,
                1068.43 ]

    def __init__(self, r: float):
        self._r = r

    def r(self, t):
        return np.full_like(t, self._r)

    def cp(self, t):
        return np.polyval(self.cp_coeff, t)

    def h(self, t):
        return np.polyval(self.h_coeff, t)
    
    def gamma(self, t):
        return self.cp(t) / (self.cp(t) - self._r)

    def phi(self, t):
        return np.polyval(self.phi_coeff[:-1], t) + self.phi_coeff[-1] * np.log(t)

    def dphi(self, t1, t2):
        return (np.polyval(self.phi_coeff[:-1], t2) - np.polyval(self.phi_coeff[:-1], t1)) + self.phi_coeff[-1] * np.log(t2 / t1)

    def pr(self, t1, t2, eff_poly):
        return np.exp(self.dphi(t1, t2) * eff_poly / self._r)

    def eff_poly(self, p1, t1, p2, t2):
        return self._r * np.log(p2 / p1) / self.dphi(t1, t2)


gpy = NumpyRealGas(287.05)

numpy_real_results = dict()
for s, n in ((1e1, 1), (1e2, 1), (1e3, 1), (1e4, 1), (1e5, 1), (2.5e5, 1), (5e5, 1), (8e5, 1), (1e6, 1), (2e6, 1), (5e6, 1), (8e6, 1), (1e7, 1)):
    merge_dicts(numpy_real_results, bench_numpy(gpy, s, repeat=10, number=n))

with open('benchs/benchs_numpy_results.json', 'w') as f:
    json.dump(numpy_real_results, f, indent=2)
