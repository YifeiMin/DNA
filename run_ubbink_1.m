
%%Set up constants for the problem
global  T kB Kb row F L;
 % nm  function of this

 
%  T=298;
%  kB=1.3*10^(-2);
%  Kb=260;
%  F=2.5;
%  L = 700;

 
 
%% Poisson Boltzman Electrostatics
moment =0; % M3
helical =0; % theta
slope = 0; % p
radius =0; % r
rho_wlc = 0;
rr=0; % dr
 
    
format long

x0_p=[0.4,10,1.5,0.01]; %x(3) is originally 2
[x_p,fval_p] = fsolve(@myfun_ubbink,x0_p);  % Call optimizer

theta_p = x_p(1);
M3_p = x_p(2);
r_p= x_p(3);
dr_p = x_p(4);

%% rho wlc from Moroz and Nelson

epsi = ( Kb*F/(kB*T)^2 - M3_p^2/(2*kB*T)^2 - 1/32)^0.5;

row1= 1-1/2*(1/epsi) + Kb*kB*T/(L*(Kb*F-M3_p^2/4)) ; % nm

row = row1;
 
slp = row*(sin(2*theta_p)/(4*pi*r_p) - (M3_p*kB*T/(8*pi*Kb*(Kb*F-(M3_p^2)/4)^(1/2))))^(-1);
%%p = row*(sin(2*theta)/(4*pi*r) - (M3*k*T/(8*pi*Kb*(Kb*F)^(1/2))))^(-1);
%%p = row*4*pi*r/(sin(2*theta));

%M3 theta r dr p (slope)
% moment = [moment;M3];
% radius = [radius;r];
% helical =[helical;theta];
% slope = [slope;p];
% rr = [rr;dr];

%% ODE Integration

%% Boltzman 

%% Analysis + Plots

%%theta_dotdot = 4*(sin(theta)).^3.*cos(theta)./(r^2) + DUpb_theta/Kb  - M3*cos(2*theta)/r/Kb;


%%Possible values of radius and angle combinations at lp

