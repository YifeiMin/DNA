
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
 
Fext1 =  0.4:0.1:1;
Fext2 =  1.1:0.1:4;

for F = Fext1
    
format long

x0= [0.1,4,1, 0.0001]; %%[0.1,1,0.5, 0.0001]; %[theta, M3, r, dr]
%%


[x,fval] = fsolve(@myfun_ubbink,x0);  % Call optimizer

theta = x(1);
M3 = x(2);
r= x(3);
dr = x(4);

%% rho wlc from Moroz and Nelson

epsi = ( Kb*F/(kB*T)^2 - M3^2/(2*kB*T)^2 - 1/32)^0.5; 

row= 1-1/2*(1/epsi) + Kb*kB*T/(L*(Kb*F-M3^2/4)) ; % nm

p = row*(sin(2*theta)/(4*pi*r) - (M3*kB*T/(8*pi*Kb*(Kb*F-(M3^2)/4)^(1/2))))^(-1); % slope
%%p = row*(sin(2*theta)/(4*pi*r) - (M3*k*T/(8*pi*Kb*(Kb*F)^(1/2))))^(-1);
%%p = row*4*pi*r/(sin(2*theta));

moment = [moment;M3];
radius = [radius;r];
helical =[helical;theta];
slope = [slope;p];
rr =[rr;dr];
end

moment1 = moment(2:length(moment));
radius1 = radius(2:length(radius));
theta1 = helical(2:length(helical));
slope1 = slope(2:length(slope));
rr1= rr(2:length(rr));

moment =0;
helical =0;
slope = 0;
radius =0;
rr=0;

for F = Fext2
    
format long

x0=[0.4,10,2,0.01];
[x,fval] = fsolve(@myfun_ubbink,x0);  % Call optimizer

theta = x(1);
M3 = x(2);
r= x(3);
dr = x(4);

%% rho wlc from Moroz and Nelson

epsi = ( Kb*F/(kB*T)^2 - M3^2/(2*kB*T)^2 - 1/32)^0.5;

row1= 1-1/2*(1/epsi) + Kb*kB*T/(L*(Kb*F-M3^2/4)) ; % nm

row = row1;
 
p = row*(sin(2*theta)/(4*pi*r) - (M3*kB*T/(8*pi*Kb*(Kb*F-(M3^2)/4)^(1/2))))^(-1);
%%p = row*(sin(2*theta)/(4*pi*r) - (M3*k*T/(8*pi*Kb*(Kb*F)^(1/2))))^(-1);
%%p = row*4*pi*r/(sin(2*theta));

moment = [moment;M3];
radius = [radius;r];
helical =[helical;theta];
slope = [slope;p];
rr = [rr;dr];

end

moment2 = moment(2:length(moment));
radius2 = radius(2:length(radius));
theta2 = helical(2:length(helical));
slope2 = slope(2:length(slope));
rr2= rr(2:length(rr));

%% ODE Integration

%% Boltzman 

%% Analysis + Plots

%%theta_dotdot = 4*(sin(theta)).^3.*cos(theta)./(r^2) + DUpb_theta/Kb  - M3*cos(2*theta)/r/Kb;


%%Possible values of radius and angle combinations at lp

