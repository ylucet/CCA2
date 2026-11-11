"""bridge_phi -- the tie residual along the sweeping family, as a polynomial in the multiplier
(port of bridgePhi.m).

For the pair g1, g2 (rescaled over d: Hessians N1, N2, linear parts m1, m2) and the point
x = xn/xd, the family M_mu = (1+mu) N1 - mu N2 has maximiser s(mu) = U(mu)/det(M_mu); Phi(mu)
is the tie value at s(mu), cleared of denominators. Returns (Phi, Dt, U1, U2), each a list of
coefficients in ASCENDING powers of mu. Exact for integer input (one path: no 2^53 wall), and
the same formula in doubles for a double x (bridgeMember's proposal filter).
"""


def _conv(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def vadd(a, b):
    n = max(len(a), len(b))
    a = list(a) + [0] * (n - len(a))
    b = list(b) + [0] * (n - len(b))
    return [x + y for x, y in zip(a, b)]


def bridgePhi(N1, m1, N2, m2, Dl, dl, kp, d, xn, xd, exact=True):
    e11 = [N1[0][0], N1[0][0] - N2[0][0]]
    e12 = [N1[0][1], N1[0][1] - N2[0][1]]
    e22 = [N1[1][1], N1[1][1] - N2[1][1]]
    Dt = vadd(_conv(e11, e22), [-x for x in _conv(e12, e12)])
    a11 = e22
    a12 = [-x for x in e12]
    a22 = e11
    w1 = [d * xn[0] - xd * m1[0], -xd * (m1[0] - m2[0])]     # xd (d x - c_mu), per component
    w2 = [d * xn[1] - xd * m1[1], -xd * (m1[1] - m2[1])]
    U1 = vadd(_conv(a11, w1), _conv(a12, w2))
    U2 = vadd(_conv(a12, w1), _conv(a22, w2))
    Phi = vadd(vadd([Dl[0][0] * x for x in _conv(U1, U1)],
                    [2 * Dl[0][1] * x for x in _conv(U1, U2)]),
               [Dl[1][1] * x for x in _conv(U2, U2)])
    Phi = vadd(Phi, [2 * xd * x for x in _conv(Dt, vadd([dl[0] * x for x in U1],
                                                        [dl[1] * x for x in U2]))])
    Phi = vadd(Phi, [2 * kp * xd * xd * x for x in _conv(Dt, Dt)])
    return Phi, Dt, U1, U2
