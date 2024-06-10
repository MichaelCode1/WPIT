"""
FLR_mod.B11
_____________________________________________________________________________________________________________________
**Description**:


_____________________________________________________________________________________________________________________
**Inputs**:

latitude: in radians

L: L shell
s: cosine of polar angle from z-axis
t: time in seconds
f0: natural frequency of field line (Hz)
Bm: wave maximum amplitude at ionoshphere (Tesla)
phi: polar angle with relation to x-axis
m: azimuthal wave number
_______________________________________________________________________________________________________________________
**Outputs**:


________________________________________________________________________________________________________________________

**Reference**:

________________________________________________________________________________________________________________________

"""

import numpy as np
from WPIT.Environment_mod import const
import WPIT.FLR_mod as flr

def B11(L, s, t, f0, Bm, phi=0, m=0):
    
    sm = flr.s0(L)
    fN = flr.fN(L)
    om0=2*np.pi*f0
    Va0=om0*L*const.Re/flr.fN(L)
    omN=Va0*flr.fN(L)/(L*const.Re)
    
    term1 = flr.h2(L, sm)*Bm*Va0**2*fN/(flr.h1(L,s)*L*const.Re*(om0**2-1/4*((omN**2-om0**2)/(4*om0))**2))
    term2 = 8*np.cos(fN*(s+s**3))*flr.f(t)
    term3 = np.cos((om0-(omN**2-om0**2)/(4*om0))*t-m*phi)*np.sinc((omN**2-om0**2)/(4*om0)*t)
    return term1*term2*term3
