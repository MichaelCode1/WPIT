"""
environment_mod.omega_plasma

**Description**:
_____________________________________________________________________________________________________________________

Routine to calculate the plasma frequency
_____________________________________________________________________________________________________________________
_____________________________________________________________________________________________________________________

**Inputs**:
_____________________________________________________________________________________________________________________

n_arg: particle number density in 𝑚−3

q_arg: particle charge in Cb

m_arg: particle mass in kg
_______________________________________________________________________________________________
_______________________________________________________________________________________________________________________

**Outputs**:
_____________________________________________________________________________________________________________________

omegap_tmp: plasma frequency in rad/s
_____________________________________________________________________________________________________________________
________________________________________________________________________________________________________________________

**Reference**:
_____________________________________________________________________________________________________________________
Parks, G. K. (1991). Physics of space plasmas. An introduction.
_____________________________________________________________________________________________________________________
________________________________________________________________________________________________________________________

"""

import numpy as np
from WPIT.Environment_mod import const

def omega_plasma(n_arg,q_arg,m_arg):
    omegap_tmp=np.sqrt((n_arg*q_arg*q_arg)/(m_arg*const.epsilon0))

    return omegap_tmp