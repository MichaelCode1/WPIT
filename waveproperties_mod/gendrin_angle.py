
import numpy as np

#####waveproperties_mod.gendrin_angle##############################

#Description:Routine to calculate the Gendrin angle
#Inputs:
# w_wave_arg: wave frequency
# wlhr_arg: lower hybrid resonance frequency
# wce_arg: electron gyrofrequency
#Outputs:
# th_gen: Gendrin angle in rad

#####################################################################

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