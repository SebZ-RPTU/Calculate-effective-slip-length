clear all

%------------------------------%
% Rectangular channel flow - Effective slip length 
%
% Description:
% Calculates the effective slip length implicitly from the analytical 
% volume flux formula for a channel flow with geometry and fluid properties 
% as defined below, assuming rectangular cross section. 
% For that, we implement a classic Newton method. 
%
% Function Inputs
%   - int_h.m
%   - vel_field_f.m
%   - volume_flux_finite_Newton_f.m
%
% Inputs:
%   - b: width of channel (in mm)
%   - h: height of channel (in mm)
%   - V: measured volume flux (in mm^3/s)
%   - dp: measured pressure drop along L (in Pa)
%   - L: length along which pressure drops (in mm)
%   - mu: dynamic viscosity of fluid (in Pa*s)
%
% Outputs:
%   - lambda_mm: effective slip length in mm
%   - lambda_um: effective slip length in microns
%
% Author:
%   Sebastian Zimmermann


%------------------------------%
%% Initialize parameters
h=0.500;
b=9.5;          

V=24.9958; 
dp=2000;             
L=10; 
mu=0.001; 




%------------------------------%
%% Initial guess for the effective slip length
lambda0=18e-3; % first guess for the effective slip length (in mm)


%------------------------------%
%% Calculate using Newton methode
NI=3;   % sets number of iterations. Higher NI lead to 
        % higher calculation times but potentially more precise results. 
hh=0.00001;
for j=1:NI;
    % value of function in step j
    F0=volume_flux_finite_Newton_f(b,h,lambda0,V,dp,mu,L);
    % now compute (estimate) of derivative;
    F0P=volume_flux_finite_Newton_f(b,h,lambda0+hh,V,dp,mu,L);
    F0M=volume_flux_finite_Newton_f(b,h,lambda0-hh,V,dp,mu,L);
    deriv=(F0P-F0M)/(2*hh);
      
    % relative error for every iteration in %
    rel_err=(abs(lambda0-round(real(lambda0-(F0/deriv)),10))/lambda0)*100 
    % iteration number
    j
    lambda0=round(real(lambda0-(F0/deriv)),10);
end;

%------------------------------%
%% Output
lambda_mm=round(lambda0,10) 
lambda_um=round(lambda0,10)*10^(3) 
    