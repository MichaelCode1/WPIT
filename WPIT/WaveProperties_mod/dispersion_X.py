"""
waveproperties_mod.dispersion_X

**Description**:
_____________________________________________________________________________________________________________________

Dispersion relation of X-mode (extra-ordinary) wave
_____________________________________________________________________________________________________________________
_____________________________________________________________________________________________________________________

**Inputs**:
_____________________________________________________________________________________________________________________

w: wave frequency

wpe: electron plasma frequency

wlh: lower hybrid resonance frequency
______________________________________________________________________________________________________________________
_______________________________________________________________________________________________________________________

**Outputs**:
_____________________________________________________________________________________________________________________

nsq_tmp: squared refractive index

kappa_tmp: wave number
______________________________________________________________________
________________________________________________________________________________________________________________________

**Reference**:
_____________________________________________________________________________________________________________________

Swanson, D. G. (2012). Plasma waves (Elsevier)
______________________________________________________________________________________________________________________
________________________________________________________________________________________________________________________

"""
from WPIT.Environment_mod import const


import numpy as np


def dispersion_X(w,wpe,wlh):
    fac1=((wpe*wpe)/(w*w))
    fac2=((w*w-wpe*wpe)/(w*w-wlh*wlh))
    nsq_tmp=1-fac1*fac2




    tmp1=w*w/(const.c_light*const.c_light)
    tmpsq=tmp1*nsq_tmp
    kappa_tmp=np.sqrt(tmpsq)
    return nsq_tmp,kappa_tmp