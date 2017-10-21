% This script computes the linking number of the rod
% z_2 and z_3 are set to be 1
%z_1 can be adjusted
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remark: all the values in this script are scaled, unless specified 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc;
clear all;
% R phi Z are centerline cylindrical coordinates
% s is the arclength with scaling
%set value of z1 first

global L  Kb Kt T b kB;
Kb = 205; %pNnm^2
Kt = 368;%pNnm^2 
T = 296.5; % Kelvin
b = Kb/Kt; %ratio of bending modulus to twisting modulus
kB = 1.3806503*10^(-2); %pNnm K^-1   Boltzmann constant


%set length of DNA molecule
bp=0.34 ;% nm % actual length of a base pair
l_bp=1000 ;% number of base pairs
L = bp*l_bp; % nm %set actual length of DNA (w/o scaling) 
lengthunit=5; % set length unit


%%%%%%%%%%%%%%%%%%%%%% set value for z1
z1= 0;
tau = sqrt((1+z1)/(1-z1));
%%

%Euler angle as a function of arclength
syms s;

varphi(s)=atan((1/(tau))*((sinh(s))/(cosh(s)))) + tau*s;
z(s)=1-(2/(1+tau^2))*(sech(s))^2;
psi(s)=atan((1/tau)*((sinh(s))/(cosh(s)))) + (3 - (2/b))*tau*s;

theta(s) = acos(z);
%%
% twist vertor as a function of arclength
kappa1(s)=diff(theta(s))*sin(psi(s))-diff(varphi(s))*sin(theta(s))*cos(psi(s));
kappa2(s)=diff(theta)*cos(psi)+diff(varphi)*sin(theta)*sin(psi);
kappa3(s)=diff(varphi)*cos(theta)+diff(psi);

kap1=@(s) eval(kappa1);
kap2=@(s) eval(kappa2);
kap3=@(s) eval(kappa3);

t = -((L/lengthunit)/2):0.01:((L/lengthunit)/2); %step size

k1=kap1(t);
k2=kap2(t);
k3=kap3(t);

plot(k1,t,'color','black'); 
hold on;
plot(k2,t,'color','b'); 
hold on;
plot(k3,t,'color','g'); 








