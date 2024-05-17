"""
environment_mod.T_bounce

**Description**:
_____________________________________________________________________________________________________________________

Calculate the bounce period of a trapped particle
_____________________________________________________________________________________________________________________
_____________________________________________________________________________________________________________________

**Inputs**:
_____________________________________________________________________________________________________________________

L_arg: L shell

v_arg: particle velocity in m/s

aeq_arg: equatorial pitch angle
_______________________________________________________________________________________________
_______________________________________________________________________________________________________________________

**Outputs**:
_____________________________________________________________________________________________________________________

bounce_tmp: Bounce period in s

_____________________________________________________________________________________________________________________
________________________________________________________________________________________________________________________

**Reference**:
_____________________________________________________________________________________________________________________
Öztürk, M. Kaan. "Trajectories of charged particles trapped in Earth’s magnetic field." American Journal of Physics 80.5 (2012): 420-428.

P. A. Sturrock, Plasma Physics: An Introduction to the Theory of Astro- physical, Geophysical, and Laboratory Plasmas (Cambridge U.P., Cam- bridge, UK, 1994)
_____________________________________________________________________________________________________________________
________________________________________________________________________________________________________________________

"""

import numpy as np
from WPIT.Environment_mod import const

def T_bounce(L_arg,v_arg,aeq_arg):
    fac1=0.117*L_arg*const.c_light/v_arg
    fac2=1-0.4635*((np.sin(aeq_arg))**(3/4))
    bounce_tmp=fac1*fac2
    return bounce_tmp
