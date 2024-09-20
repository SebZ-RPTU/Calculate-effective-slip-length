function p=int_h(x,b,h,lambda_eff,dp,mu,L)

%------------------------------%
% Description:
% Calculates the y-integral for every given x for the velocity field.
%
% Function Inputs
%   - vel_field_f

%------------------------------%
%% integral
N1=1000; % the higher this number, the finer the mesh for the integral.
hh = h/N1;

sumh=0;
yy=0;
sumh = sumh+(hh/2)*vel_field_f(x+i*yy,b,h,lambda_eff,dp,mu,L);
    for k=1:N1-1; 
        yy=0+k*hh;
        sumh=sumh+hh*vel_field_f(x+i*yy,b,h,lambda_eff,dp,mu,L); 
    end;

% y=h doesn't contribute to the integral. Therefore no last part of the
% sum. 

%------------------------------%
%% function value
p = sumh;