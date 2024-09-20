function p=volume_flux_finite_Newton_f(b,h,lambda_eff,V,dp,mu,L)

%------------------------------%
% Description:
% Calculates the volume flux for the rectangular channel. For that, we
% numercially integrate the velocity field w(x,y) for the whole 
% cross- section.  
%
% Function Inputs
%   - int_h.m

%------------------------------%
%% Initialize integration
% setting mesh resolution
Nstep=1000; % the higher this number, the finer the mesh for the integral. 
step = b/(2*Nstep);

% setting step size
XX = 0:step:(b/2);
bb = b/(2*length(XX));

%------------------------------%
%% Calculate integral (using trapeziodal method) 
% We wrote a function that calculates the y-integral for every 
% given x (see int_h.m).
% Because of that, here, we only have to add up all the x-values 
% in a for-loop.
sum=0;
sum = 0+(bb/4)*int_h(0,b,h,lambda_eff,dp,mu,L); 
% only take half of the contribution (hence 'bb/4') of x=0, since here, we only calculate
% for the domain x>=0 (due to x-axis symmetry)
for j=1:(length(XX)-1);
    sum = sum+bb*int_h(XX(j),b,h,lambda_eff,dp,mu,L);
end;
sum = sum+0; % integral must be zero along the wall. 


%------------------------------%
%% Function value
% this is the analytical result for the volume flux V!
p = (2*sum)-V;