"""
FLR_mod.SnTemp
_____________________________________________________________________________________________________________________
**Description**:


_____________________________________________________________________________________________________________________
**Inputs**:

L: L shell
s: cosine of polar angle from z-axis
A: 
_______________________________________________________________________________________________________________________
**Outputs**:


________________________________________________________________________________________________________________________

**Reference**:

________________________________________________________________________________________________________________________

"""

import numpy as np
from WPIT.Environment_mod import const
import WPIT.FLR_mod as flr

def SnTemp(L, s):

    sm = flr.s0(L)
    fN = flr.fN(L)
    A = L**2*const.Re**2/(sm+sm**3)**(1/2)
    
    return A*np.sin(fN*(s+s**3))
