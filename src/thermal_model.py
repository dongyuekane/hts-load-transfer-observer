"""Synthetic two-temperature cable model and the balanced two-bus load model.

The equivalent-strip heating calculation and power-flow equations are from
this study's synthetic simulation. Electrical and thermal coefficients are
simulation assumptions.
"""

from functools import lru_cache
import math

import numpy as np

CC=180000.; CN=1020000.; H=800.; G=750.; Q0=3000.
QR=6000.; IC77=6000.; TS=70.; TC=90.; IR=2200.
def strip(v):
    v=np.minimum(np.maximum(v,0.),.999999)
    z=(1-v)*np.log1p(-v)+(1+v)*np.log1p(v)-v*v
    return np.where(v<.01,v**4/6+v**6/15+v**8/28,z)
NORM=float(strip(np.sqrt(2)*IR/IC77))
def ic(T):return IC77*np.maximum((TC-T)/13,.001)
def heat(I,T):
    ci=ic(T)
    return QR*(ci/IC77)**2*strip(np.sqrt(2)*I/ci)/NORM
def bound(I):return np.minimum(85.,90.-13*np.sqrt(2)*I/(.8*IC77))

def powerflow(p_mw, branch):
    """Exact balanced constant-PQ two-bus AC solution, 100 MVA / 35 kV base.
    Each trunk supplies isolated aggregated load blocks; no source paralleling.
    Current returned is sending-end RMS line current. No shunt admittance.
    """
    z = .001+ .005j if branch == 'hts' else .025 + .06j
    s = p_mw / 100 * (1+1j*math.tan(math.acos(.95)))
    v = 1+0j
    for _ in range(100):
        cur = np.conj(s/v)
        vn = 1-z*cur
        if abs(vn-v)<1e-12:
            v=vn
            break
        v=vn
    cur = np.conj(s/v)
    # Independent branch KCL residual, not the identity V * conj(S*/V*) = S.
    mismatch = abs((1-v)/z-cur)
    return abs(v), abs(cur)*100e6/(math.sqrt(3)*35e3), (abs(cur)**2*z.real*100), mismatch


@lru_cache(None)
def electrical(load):
    p=load/4
    h=[powerflow(n*p,'hts') for n in range(5)]
    b=[powerflow(n*p,'backup') for n in range(5)]
    return np.array([v[1] for v in h]),np.array([v[0] for v in b])

def rhs(state, current, cooling, capacity_lower, capacity_upper, loss=1., conductance=1.):
    conductor, coolant = state[..., 0], state[..., 1]
    flux = H * conductance * (conductor - coolant)
    flow = np.stack(
        [heat(current, conductor) * loss - flux, flux + Q0 - cooling * G * (coolant - TS)],
        axis=-1,
    )
    capacities = np.array([CC, CN]) * np.where(flow >= 0, capacity_lower, capacity_upper)
    return flow / capacities
