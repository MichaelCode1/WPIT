"""
FLR_mod.B2
_____________________________________________________________________________________________________________________
**Description**:


_____________________________________________________________________________________________________________________
**Inputs**:

L: L shell
s: cosine of polar angle from z-axis
t: time in seconds
f0: natural frequency of field line (Hz)
Bm: wave maximum amplitude at ionoshphere (Tesla)
phi: polar angle with relation to x-axis
m: wave number
_______________________________________________________________________________________________________________________
**Outputs**:


________________________________________________________________________________________________________________________

**Reference**:

________________________________________________________________________________________________________________________

"""

import numpy as np
from WPIT.Environment_mod import const
import WPIT.FLR_mod as flr

def B2(L, s, t, f0, Bm, phi=0, m=0):
    
    sm = flr.s0(L)
    fN=flr.fN(L)
    f=flr.f(t)
    om0=2*np.pi*f0
    Va0=om0*L*const.Re/fN
    omN=Va0*fN/(L*const.Re)
    
    term1 = flr.h2(L, sm)/flr.h2(L,s)*Bm*np.sin(fN*(s+s**3))*flr.f(t)
    term2 = np.cos((om0+(omN**2-om0**2)/(4*om0))*t - m*phi)
    term3 = np.sinc((omN**2-om0**2)/(4*om0)*t)
    
    return term1*term2*term3
