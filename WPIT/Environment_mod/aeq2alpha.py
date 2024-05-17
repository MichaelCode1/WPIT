"""
environment_mod.aeq2alpha

**Description**:
_____________________________________________________________________________________________________________________

Routine to translate equatorial pitch angle to local pitch angle

_____________________________________________________________________________________________________________________
_____________________________________________________________________________________________________________________

**Inputs**:
_____________________________________________________________________________________________________________________

L_arg: L shell

lambda_arg: magnetic latitude in rad

aeq_arg: equatorial pitch angle in rad

_____________________________________________________________________________________________________________________
_______________________________________________________________________________________________________________________

**Outputs**:
_____________________________________________________________________________________________________________________

alpha0: local pitch angle in rad

_____________________________________________________________________________________________________________________
________________________________________________________________________________________________________________________

**Reference**:

_____________________________________________________________________________________________________________________

_____________________________________________________________________________________________________________________
________________________________________________________________________________________________________________________

"""

import numpy as np
from WPIT.Environment_mod.Bmag_dipole import Bmag_dipole


def aeq2alpha(L_arg,lambda_arg,aeq_arg):
    Blam0=Bmag_dipole(L_arg,lambda_arg)
    Beq0=Bmag_dipole(L_arg,0)
    salpha0=np.sin(aeq_arg)*np.sqrt(Blam0/Beq0)
    alpha0=np.arcsin(salpha0)
    
    return alpha0
