from pythermo import PolyGas8, mix8
import copy

class AirFarWar(PolyGas8):
    """Gas model of air with possible water and/or fuel."""

    def __init__(self, far, war):
        cps = [
            (
                2.35822e-20,
                -1.79613e-16,
                4.70972e-13,
                -3.3094e-10,
                -6.27984e-07,
                0.00123785,
                -0.521742,
                1068.43,
            ),
            (0.0, 0.0, 0.0, 0.0, 0.0, -0.000527394, 2.3021158, 1285.0007),
            (0.0, 0.0, -3.86134e-14, 3.345e-10, -1.14179e-6, 0.00173478, -0.484304, 1884.57),
        ]
        h0s = [0.0, -42785.8, 0.0]
        rs = [287.05287, 287.05287, 461.522]
        w = [1.0, far, war]
        mix_far_war = mix8(cps, h0s, rs, w)

        super().__init__(mix_far_war)


g = AirFarWar(0., 0.)
print(g.h(288))
print(g.h(300))

g2 = copy.deepcopy(g)
print(g2.h(288))
print(g2.h(300))