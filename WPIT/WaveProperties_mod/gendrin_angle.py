"""
waveproperties_mod.gendrin_angle

**Description**:
_____________________________________________________________________________________________________________________

Routine to calculate the Gendrin angle
_____________________________________________________________________________________________________________________
_____________________________________________________________________________________________________________________

**Inputs**:
_____________________________________________________________________________________________________________________

w_wave_arg: wave frequency

wlhr_arg: lower hybrid resonance frequency

wce_arg: electron gyrofrequency
______________________________________________________________________________________________________________________
_______________________________________________________________________________________________________________________

**Outputs**:
_____________________________________________________________________________________________________________________

th_gen: Gendrin angle in rad

______________________________________________________________________
________________________________________________________________________________________________________________________

**Reference**:
_____________________________________________________________________________________________________________________

Bortnik, J., Inan, U. S., and Bell, T. F. (2006). Landau damping and resultant unidirectional propagation of
chorus waves. Geophysical research letters 33
______________________________________________________________________________________________________________________
________________________________________________________________________________________________________________________

"""
import numpy as np



def gendrin_angle(w_wave_arg,wlhr_arg,wce_arg):
    if w_wave_arg<wlhr_arg:
        th_gen=np.nan
    if w_wave_arg==wlhr_arg:
        th_gen=np.deg2rad(90)
    if wlhr_arg<w_wave_arg<0.5*wce_arg:
        th_gen=np.arccos(2*(w_wave_arg*w_wave_arg-wlhr_arg*wlhr_arg)/(w_wave_arg*wce_arg))
    if 0.5*wce_arg<w_wave_arg<wce_arg:
        th_gen=0

    return th_gen