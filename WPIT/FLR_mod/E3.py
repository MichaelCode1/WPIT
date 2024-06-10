"""
FLR_mod.E3
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

def E3(L, s, t, f0, Bm, phi=0, m=0):
   
    sm = (1 - 1/L)**(1/2)
    om0=2*np.pi*f0
    Va0=om0*L*const.Re/flr.fN(L)
    omN=Va0*flr.fN(L)/(L*const.Re)
    
    term1 = flr.h2(L,sm)*Bm*Va0**2*flr.fN(L)/(flr.h3(L,s)*L**2*const.Re**2*(om0+(omN**2-om0**2)/(4*om0**2)))
    term2 = np.cos(flr.fN(L)*(s+s**3)) * flr.f(t)
    term3 = np.sin((om0+(omN**2-om0**2)/(4*om0))*t-m*phi)*np.sinc((omN**2-om0**2)/(4*om0)*t)
    
    return term1*term2*term3