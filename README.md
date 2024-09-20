# Calculate-effective-slip-length
A tool to calculate the effective slip length for laminar flow in rectangular channels. 

Description:
Calculates the effective slip length implicitly from the analytical 
volume flux formula for a channel flow with geometry and fluid properties 
as defined below, assuming rectangular cross section. 
For that, we implement a classic Newton method. 

Function Inputs
  - int_h.m
  - vel_field_f.m
  - volume_flux_finite_Newton_f.m

 Inputs:
  - b: width of channel (in mm)
  - h: height of channel (in mm)
  - V: measured volume flux (in mm^3/s)
  - dp: measured pressure drop along L (in Pa)
  - L: length along which pressure drops (in mm)
  - mu: dynamic viscosity of fluid (in Pa*s)

 Outputs:
  - lambda_mm: effective slip length in mm
  - lambda_um: effective slip length in microns

Author: Sebastian Zimmermann

Original publication: [...]
