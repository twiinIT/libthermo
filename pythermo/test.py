from pythermo import IdealGas
from pythermo import PolyGas8
from pythermo import PolyGas8, PolyIdealGasProps8
import numpy as np

gas = IdealGas(288.058, 1004.4)
print(gas.gamma(288.15))

props = PolyIdealGasProps8(
    [
        2.35822e-20,
        -1.79613e-16,
        4.70972e-13,
        -3.3094e-10,
        -6.27984e-07,
        0.00123785,
        -0.521742,
        1068.43,
    ],
    0,
    287.058,
)
gas = PolyGas8(props)
print(gas.gamma(288.15))
print(
    gas.pressure_ratio(
        np.full(1_000_000_000, 288.15),
        np.full(1_000_000_000, 400.15),
        np.full(1_000_000_000, 0.9),
    ).shape
)
