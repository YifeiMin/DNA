% This script computes the linking number of the rod
% z_2 and z_3 are set to be 1
%z_1 can be adjusted
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remark: all the values in this script are scaled (including curvature twisting etc.), unless specified 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear all;
% R phi Z are centerline cylindrical coordinates
% s is the arclength with scaling
%set value of z1 first

global massDNA Kb Kt T b kB F ;
F = 2.5 ; % pN
Kb = 205; %pNnm^2
Kt = 368;%pNnm^2 
T = 296.5; % Kelvin
b = Kt/Kb; %ratio of twisting modulus to bending modulus
kB = 1.3806503*10^(-2); %pNnm K^-1   Boltzmann constant

%set length and mass of DNA molecule
bp=0.34 ;% nm % actual length of a base pair
l_bp=1900 ;% number of base pairs
L = bp*l_bp; % nm %set actual length of DNA (w/o scaling) 
%massDNA = l_bp*650*1.67*(10^(-24)) ; % gram   650 dalton is the mass for one base pair 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global z1; 
z1= 0.7; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%set unit for mass, length, and time
%Munit = rho*sqrt(A*I1); % gram  Munit=[M] is the scaled unit of mass 

Lunit= sqrt((2*Kb)/(F*(1-z1))); % nm  length unit = [L]
%Tunit= sqrt((Munit/Kb)*(Lunit^3)); % s  time unit = [T]   
Torqueunit = sqrt(Kb*F);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% set value for tau
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

t = -((L/Lunit)/2):0.01:((L/Lunit)/2); %step size

k1=kap1(t);
k2=kap2(t);
k3=kap3(t);

%plot twist vector
figure;
plot(k1,t,'color','black','linewidth',1); 
hold on;
plot(k2,t,'color','b','linewidth',1); 
hold on;
plot(k3,t,'color','g','linewidth',1); 
%%

fR =(2/(1+(tau)^2))*sech(s) ;

fphi = tau*s - pi/2;

%transform into cartisian coordinates
fX = fR*cos(fphi);
fY = fR*sin(fphi);
fZ = s- (2 / (1+tau^2))*(sinh(s)/cosh(s)) ;


Xs = @(s) eval(fX);
Ys = @(s) eval(fY);
Zs = @(s) eval(fZ);

%plot centerline 

X=Xs(t);
Y=Ys(t);
Z=Zs(t);

figure;
plot(t,Zs(t));

figure;
plot3(X,Y,Z, 'linewidth',2);
axis([-5 5 -5 5 -20 20 ])

M3actual = sqrt((1+z1)*(1+1)*(Kb*F)); %pNnm
