%% This script: Plots Helical Angle and Radius versus Salt Concentration

% Salt Concentration is 1[M]. Deybe length and v is changed accordingly.
%compute straight state and inter state only
clc;
clear all;
% R phi Z are centerline cylindrical coordinates
% s is the arclength with scaling
%set value of z1 first

global L Kb Kt T b kB F bp l_bp z1 lB cr Wr dn a kd v eta nlp;
F = 2 ; % pN
Kb = 205; %pNnm^2 % 205 is also a chioce
Kt = 389;%pNnm^2 
T = 296.5; % Kelvin
b = Kt/Kb; %ratio of twisting modulus to bending modulus
kB = 1.3806503*10^(-2); %pNnm K^-1   Boltzmann constant

lB= 0.715; %nm
cr = 1/(2^(8/3)); % parameter control entropy
eta = cr/(pi)^(2/3);
Wr = 1; % parameter control size loop
dn = 0; % delta n at the transiton 0.7 for short and 1 for long 
a=1;  % parameter control fluctuations
y= 1;


%set length and mass of DNA molecule
bp=0.34 ;% nm % actual length of a base pair
l_bp=6000 ;% number of base pairs
L = bp*l_bp; % nm %set actual length of DNA (w/o scaling) 

nlp=11;

z1range = -0.999:0.01:-0.8;
V_3vari = [0];
V_strvari=[0];
lknum=[0];

vecr=[0]; %helical raduis
vectheta=[0]; %helical angle
ic=[0];

for ioncon = 0.1 : 0.05 : 0.8
v = 2.46+2.38*11*ioncon ; %0.72*8.231260437 ;  % 0.72 47.7 ; %8.06 ;%% 1.97; % nm^-1
kd =1/(0.305/sqrt(ioncon));%1/0.351; %1/0.8; %%1/3.07; %mM %reciprocal of lambda_D, eqn(18)
%%2.2 kbp variables
run_ubbink_1  %%   (Amoment)co = 150 mM  Marko and Siggia Electrostatics and Entropy 
theta_p_deg = theta_p*180/pi;

ic=[ic,ioncon];
vecr=[vecr,r_p];
vectheta=[vectheta, theta_p_deg]

end

figure
plot(ic(2:length(ic)) , vecr(2:length(vecr)), 'b','LineWidth',2);

figure
plot(ic(2:length(ic)) , vectheta(2:length(vectheta)), 'r','LineWidth',2);



%% lk num and energy as a function of Lp under plectonemic regime
A = (Kb/kB/T)^(1/3);
entropy_p = cr*kB*T/(A*(dr_p^(2/3))) + eta*kB*T/(A*((r_p*cot(theta_p)))^(2/3));
C = (1/2)*kB*T*v^2*lB*(pi/kd/r_p)^(1/2);
g = 1 + 0.83*(tan(theta_p))^2+0.86*(tan(theta_p))^4;
Upb_p = C*g*exp(2*kd^2*dr_p^2-2*kd*r_p) + entropy_p ; 
K_p = ((Kb*F-((M3_p)^2)/4)^(1/2))/(kB*T);  % K in eqn(6)
Gflu_p = ((kB*T)^2*(K_p)/Kb)*(1-1/(4*(K_p))-1/(64*(K_p)^2));
Lkhat_p = M3_p^2/(4*Kb*K_p); 
Lo_p = 4*(Kb/F)^(1/2);
kappa_p = ((sin(theta_p))^2)/r_p;





% 
% V_p0=[0];
% n_p0=[0];
% 
% for Lp_p= (0.03*L):(0.001*L): (0.5*L)
%     
% V_p1 = L*(- F + (M3_p)^2/(2*Kt) + Gflu_p + Lkhat_p ) + Lo_p*( 2*F + a*( - Gflu_p - Lkhat_p) ) + (Lp_p)*(F + Upb_p - Gflu_p - Lkhat_p )+  1/2*(kappa_p^2*Kb*Lp_p) ; %energy of plectonemic regime
% 
% n_p1 = M3_p*L/(2*pi)*(1/Kt + 1/(4*Kb*K_p)) - (a*Lo_p+Lp_p)*M3_p/(2*pi)*(1/(4*Kb*K_p)) + (Wr-dn) + sin(2*theta_p)/(4*pi*r_p)*Lp_p; %link of plectoneme
% 
% V_p0=[V_p0 , V_p1 ];
% n_p0=[n_p0 , n_p1 ];
% 
% end
% 
% 
% 
% figure
% plot(lknum(2:length(lknum)) , V_3vari(2:length(V_3vari)), 'b','LineWidth',2);
% hold on;
% plot(lknum(2:length(lknum)) , V_strvari(2:length(V_strvari)), 'r','LineWidth',2);
% hold on;
% plot(n_p0(2:length(n_p0)), V_p0(2:length(V_p0)),'g', 'linewidth' ,2);
% xlabel('Linking num','fontsize',14,'fontweight','b'); ylabel('energy[pNnm]','fontsize',14,'fontweight','b');





