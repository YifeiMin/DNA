clear all;
F=3.5;
M3 = 19.4;
r=1.3017;
theta=0.3697;
Kb = 205;
kB= 1.3806503*10^(-2);
T = 296.5;

K=sqrt(Kb*F-((M3)^2)/4)/(kB*T);

dndLp = sin(2*theta)/(4*pi*r)-M3/(8*pi*Kb*K); 
dLp = 4*pi*r/sin(theta);
dn=dndLp*dLp;