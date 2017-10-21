%% final script for : find transition point 1,2, then plot stair step graph
%% F and ion concentration are fixed in this script

clc;
clear all;
% R phi Z are centerline cylindrical coordinates
% s is the arclength with scaling
%set value of z1 first

%% set up constants
global L Kb Kt T b kB F bp l_bp lB cr Wr dn a kd v eta ;
F = 3.5 ; % pN
concen=1; % molar, ion concentration 
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
v = 2.46+2.38*10^(-2)*(concen*10^3); %0.72*8.231260437 % 0.72 47.7 ; %8.06 ;%% 1.97; % nm^-1
kd =sqrt(concen)/0.305;%1/0.351; %1/0.8; %%1/3.07; %mM %reciprocal of lambda_D, eqn(18)

%set length and mass of DNA molecule
bp=0.34 ;% nm % actual length of a base pair
l_bp=5882 ;% number of base pairs
L = bp*l_bp; % nm %set actual length of DNA (w/o scaling) 


%% find critical link and moment at first transition
run_transition1; %% find critical link and moment at first transition 
% Variable Output: Mcric_s  Lp0_s K_s ncric_s

%% plot M vs n curve BEFORE transition1
M3vec_s=[0];
Lkvec_s=[0];

for M3 = 0:0.1:Mcric_s
    M3vec_s=[M3vec_s, M3];
    K=sqrt(Kb*F - M3^2/4)/(kB*T);
    
    Lk=M3/(2*pi)*((1/Kt)+(1/(4*Kb*K)))*L;
    Lkvec_s=[Lkvec_s,Lk];
end

%plot M vs n before transition
% figure
% plot( Lkvec_s(2:length(Lkvec_s)),M3vec_s(2:length(M3vec_s)),'--','LineWidth',1); %2
% xlabel('n (turn)');
% ylabel('M3');
% title('M3 vs n in before transition');

approximate_Mslope_s = ( M3vec_s(length(M3vec_s)-10) -  M3vec_s(length(M3vec_s)-11) )/(Lkvec_s(length(Lkvec_s)-10)-Lkvec_s(length(Lkvec_s)-11));   

%%
run_ubbink_1;
% variable output:   theta_p  M3_p  r_p  dr_p 

K_p = sqrt(Kb*F - M3_p^2/4)/(kB*T);
Mmax1=sqrt(187*sqrt(F)*2/L/((1/Kt)+(1/(4*Kb*K_p))+(M3_p^2)/(16*Kb*K_p^3*(kB^2)*(T^2)))   +M3_p^2); 

diffn=(Mmax1-M3_p)*L/(2*pi)*((1/Kt)+(1/(4*Kb*K_p))+M3_p^2/(16*Kb*K_p^3*kB^2*T^2));

figure
plot( Lkvec_s(2:length(Lkvec_s)),M3vec_s(2:length(M3vec_s)),'LineWidth',1); %2
xlabel('n (turn)');
ylabel('Torque(pNnm)');
title(['F=',num2str(F),'pN , c_0=1molar']);

axis([0 70 0 35])
hold on
for n_dum = ncric_s: diffn: (ncric_s+18*diffn)
   pt1=[n_dum,n_dum+diffn];
   pt2=[M3_p,Mmax1];
   
   plot(pt1,pt2,'b','LineWidth',1)
   hold on
   
end


