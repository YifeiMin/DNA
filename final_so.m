%compute straight state and inter state only
clc;
clear all;
% R phi Z are centerline cylindrical coordinates
% s is the arclength with scaling
%set value of z1 first

global L Kb Kt T b kB F bp l_bp z1;
F = 2.5 ; % pN
Kb = 205; %pNnm^2 % 205 is also a chioce
Kt = 389;%pNnm^2 
T = 296.5; % Kelvin
b = Kt/Kb; %ratio of twisting modulus to bending modulus
kB = 1.3806503*10^(-2); %pNnm K^-1   Boltzmann constant

%set length and mass of DNA molecule
bp=0.34 ;% nm % actual length of a base pair
l_bp=2200 ;% number of base pairs
L = bp*l_bp; % nm %set actual length of DNA (w/o scaling) 


z1range = -0.999:0.01:-0.8;
V_3vari = [0];
V_strvari=[0];
lknum=[0];

for z1 = z1range

energy_inter_z1 ;
energy_straight ;
V_3vari=[ V_3vari , V_3] ;
V_strvari=[V_strvari; V_straight] ;
lknum=[lknum; Lk] ;

end

figure
plot(z1range,V_strvari(2:length(V_strvari)),'r','LineWidth',2);
hold on;
plot(z1range,V_3vari(2:length(V_3vari)), 'b','LineWidth',2);
xlabel('z1','fontsize',14,'fontweight','b'); ylabel('energy[pNnm]','fontsize',14,'fontweight','b');


figure
plot(lknum(2:length(lknum)) , V_3vari(2:length(V_3vari)), 'b','LineWidth',2);
hold on;
plot(lknum(2:length(lknum)) , V_strvari(2:length(V_strvari)), 'r','LineWidth',2);
xlabel('Linking num','fontsize',14,'fontweight','b'); ylabel('energy[pNnm]','fontsize',14,'fontweight','b');
