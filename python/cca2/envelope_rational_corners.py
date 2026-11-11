"""envelope_rational_corners -- move every CORNER face whose dual point is RATIONAL from the
bridge list to the quadratic list, where it is an ordinary affine face
(port of envelopeRationalCorners.m).
"""
from math import gcd, isqrt

from . import rat_q


def envelopeRationalCorners(env):
    env = {'quad': list(env['quad']), 'bridge': list(env['bridge'])}
    keep = [True] * len(env['bridge'])
    for k, c in enumerate(env['bridge']):
        if c.get('kind') != 'corner' or not c.get('con') or 'sigma' not in c:
            continue
        sn, sd, ok = _rationalCorner(c['sigma'])
        if not ok:
            continue
        n = c['src']['n']
        d = c['src']['d']
        N = [[n[4], n[5]], [n[5], n[6]]]
        m = [n[7], n[8]]
        qq = n[9]
        g0 = rat_q.quad(sn, N, sn) + 2 * sd * rat_q.dot(m, sn) + 2 * sd ** 2 * qq
        num = [0, 0, 0, 0, 0, 0, 0, 2 * d * sd * sn[0], 2 * d * sd * sn[1], -g0]
        den = 2 * d * sd ** 2
        gg = gcd(den, gcd(gcd(num[7], num[8]), num[9]))
        env['quad'].append({'num': [x // gg for x in num], 'den': den // gg, 'con': c['con'],
                            'pair': c['pairIndex']})
        keep[k] = False
    env['bridge'] = [b for b, kp in zip(env['bridge'], keep) if kp]
    return env


def _rationalCorner(sg):
    al, be, ga = sg['quad'][0], sg['quad'][1], sg['quad'][2]
    if sg['root'] == 0 or al == 0:
        if be == 0:
            return [0, 0], 1, False
        p = -ga
        q = be                                      # sigma = -ga/be
    else:
        D = be * be - 4 * al * ga
        if D < 0:
            return [0, 0], 1, False
        r = isqrt(D)
        if r * r != D:
            return [0, 0], 1, False                 # an irrational corner stays a bridge
        p = -be + sg['root'] * r
        q = 2 * al
    sn = [sg['s0n'][0] * q + p * sg['w'][0], sg['s0n'][1] * q + p * sg['w'][1]]
    sd = sg['na'] * q
    if sd < 0:
        sn = [-sn[0], -sn[1]]
        sd = -sd
    gs = gcd(gcd(sn[0], sn[1]), sd)
    return [sn[0] // gs, sn[1] // gs], sd // gs, True
