# GravityWave2D_Propagation
2D numerical simulation of atmospheric gravity wave propagation using FDTD method.

Fully nonlinear, compressible 2D fluid equations are solved with gravity and simple molecular viscosity & thermal diffusivity terms.
The main (advective) equations are solved using dimensionally-split Lax-Wendroff, and the diffusive equations are solved using an explicit method in a time-split manner (Godunov splitting).

This general 2D propagation model can be used to link wave sources and GW effects in the upper atmosphere, and can be applied to other planets. Both acoustic and gravity waves are captured in this model's physics. 
