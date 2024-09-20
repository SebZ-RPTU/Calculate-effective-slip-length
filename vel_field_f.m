function p=vel_field_f(z,b,h,lambda_eff,dp,mu,L)

%------------------------------%
% Description:
% Calculates the velocity field for the rectangular channel flow
% analytically, as defined in the paper. This will later be integrated to
% receive the corresponding volume flux. 
%
% Function Inputs
%   - none

%------------------------------%
%% Initialize parameter for infinite sums
N1=19; % as the function contains infinite sums, we 
       % need to truncate at a certain term. Lower N1 benefits 
       % stability if b>>h. 
x=real(z);
y=imag(z);

%------------------------------%
%% Calculate the constant Omega
% numerator
sum_omega_num=0;
for k=1:2:N1;
    sum_omega_num=sum_omega_num+(1/k^3)*(1-(cosh((k*pi*x)/(2*h))/cosh((k*pi*b)/(4*h))))*sin((k*pi)/2);
end;
% demoninator
sum_omega_denom=0;
for k=1:2:N1;
    sum_omega_denom=sum_omega_denom+(1/k^2)*(1-(cosh((k*pi*x)/(h))/cosh((k*pi*b)/(2*h))));
end;
% finally get Omega in terms of lambda_eff
Omega = (1/lambda_eff)*(sum_omega_num/sum_omega_denom);


%------------------------------%
%% Calculate velocity field w
% sum1
sum1=0;
for k=1:2:N1;
    sum1=sum1+(1/k^3)*(1-(cosh((k*pi*x)/(h))/cosh((k*pi*b)/(2*h))))*sin((k*pi*y)/h);
end;
% sum2
sum2=0;
for k=1:2:N1;
    sum2=sum2+(1/k^3)*(1-(cosh((k*pi*x)/(2*h))/cosh((k*pi*b)/(4*h))))*sin((k*pi)/2)*cos((k*pi*y)/(2*h));
end;

%------------------------------%
%% Function return value
% this is the velocity field w(x,y)!
p=(dp/(mu*L))*((16*h^2)/pi^3)*(  ((Omega*h)/(4*Omega*h+pi))*sum1  +  (pi/(4*Omega*h+pi))*sum2   );


