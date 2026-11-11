"""_poly_sign_nb -- the Bernstein sign test of poly_sign_on, compiled with numba.

Same arithmetic as poly_sign_on._bernstein (the change of variable x = x0 + h u on the
coefficients, then the power-to-Bernstein matrix, on both sides) and the same margin
1e-10 * max(max|B|, eps). Only the order of the floating-point sums differs from numpy's
matrix products, i.e. rounding; every decision is taken behind the 1e-10 relative margin.
"""
import numpy as np
from numba import njit


@njit(cache=True)
def _shift(n, x0, h, Cb):
    S = np.zeros((n, n))
    hk = 1.0
    for k in range(n):
        for i in range(k, n):
            S[k, i] = Cb[i, k] * x0 ** (i - k) * hk
        hk *= h
    return S


@njit(cache=True)
def bernSign(P, x0, x1, y0, y1, Cb, M1, M2):
    n1, n2 = P.shape
    Sx = _shift(n1, x0, x1 - x0, Cb)
    Sy = _shift(n2, y0, y1 - y0, Cb)
    # A = Sx P Sy^T, B = M1 A M2^T
    T = np.zeros((n1, n2))
    for i in range(n1):
        for j in range(n2):
            acc = 0.0
            for k in range(n1):
                acc += Sx[i, k] * P[k, j]
            T[i, j] = acc
    A = np.zeros((n1, n2))
    for i in range(n1):
        for j in range(n2):
            acc = 0.0
            for k in range(n2):
                acc += T[i, k] * Sy[j, k]
            A[i, j] = acc
    U = np.zeros((n1, n2))
    for i in range(n1):
        for j in range(n2):
            acc = 0.0
            for k in range(n1):
                acc += M1[i, k] * A[k, j]
            U[i, j] = acc
    m = 0.0
    B = np.zeros((n1, n2))
    for i in range(n1):
        for j in range(n2):
            acc = 0.0
            for k in range(n2):
                acc += U[i, k] * M2[j, k]
            B[i, j] = acc
            if abs(acc) > m:
                m = abs(acc)
    tol = 1e-10 * max(m, 2.0 ** -52)
    pos = True
    neg = True
    for i in range(n1):
        for j in range(n2):
            if not B[i, j] > tol:
                pos = False
            if not B[i, j] < -tol:
                neg = False
    if pos:
        return 1
    if neg:
        return -1
    return 0
