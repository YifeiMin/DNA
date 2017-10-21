
clc;
clear all;
F_ep = 3.5 ; % pN
M_ep = 22.7 ; %pNnm max moment (stair-step) from experiment
omg1 = 2 ; % rad/sec 
fc1=1 ; % sec^-1
omg2 = 1 ; % rad/sec 
fc2=1 ; % sec^-1
omg3 = 1 ; % rad/sec 
fc3=5 ; % sec^-1
omg4 = 1 ; % rad/sec 
fc4=10 ; % sec^-1

global L Kb Kt T b kB F bp l_bp lB cr Wr dn a kd v eta ;
Kb = 205; % pNnm^2
Kt = 389;%pNnm^2 
b = Kt/Kb; %ratio of twisting modulus to bending modulus
lB= 0.715; %nm
cr = 1/(2^(8/3)); % parameter control entropy
eta = cr/(pi)^(2/3);
Wr = 1; % parameter control size loop
dn = 0; % delta n at the transiton 0.7 for short and 1 for long 
a=1;  % parameter control fluctuations
y= 1;
v = 2.46+2.38*10^(-2)*(10^3) ;  % 0.72 47.7 ; %8.06 ;%% 1.97; % nm^-1
kd =sqrt(1)/0.305;%1/0.351; %1/0.8; %%1/3.07; %mM %reciprocal of lambda_D, eqn(18)

bp=0.34 ;% nm % actual length of a base pair
l_bp=5882 ;% number of base pairs
L = bp*l_bp; % nm %set actual length of DNA (w/o scaling) 

kB = 1.3806503*10^(-2); %pNnm K^-1   Boltzmann constant
T = 296.5; % Kelvin
c1 = (M_ep - (kB*T)/(2*pi)*log(omg1/(4*pi*fc1)) )/sqrt(F_ep); %omg/fc = 2
c2 = (M_ep - (kB*T)/(2*pi)*log(omg2/(4*pi*fc2)) )/sqrt(F_ep); %1
c3 = (M_ep - (kB*T)/(2*pi)*log(omg3/(4*pi*fc3)) )/sqrt(F_ep); %1/5
c4 = (M_ep - (kB*T)/(2*pi)*log(omg4/(4*pi*fc4)) )/sqrt(F_ep); %1/10

F=3.5;
format long 
M_0 = 25; % initial search value pNnm
[Mcrit,fval] = fsolve(@findMcrit1,M_0); % Mcrit is just the critical moment before transition



Fvec1=[0];
Mvec1=[0];
Fvec2=[0];
Mvec2=[0];
Fvec3=[0];
Mvec3=[0];
Fvec4=[0];
Mvec4=[0];


Mpvec1=[0]; 
Mmaxvec1=[0];
dnvec=[0];
diffMvec=[0];
slpvec=[0];

for Frange = 1:0.1:10
    F = Frange; 
    run_ubbink_1; 
    Mpvec1=[Mpvec1,M3_p];
    
    K=sqrt(Kb*F-M3_p^2/4)/(kB*T);
    
    
    Mmax1=sqrt(187*sqrt(Frange)*2/L/((1/Kt)+(1/(4*Kb*K))+(M3_p^2)/(16*Kb*K^3*(kB^2)*(T^2)))   +M3_p^2); % critical moment 2 (stair step)
    diffn=(Mmax1-M3_p)*L/(2*pi)*((1/Kt)+(1/(4*Kb*K))+M3_p^2/(16*Kb*K^3*kB^2*T^2)); % see report_Mmax_vs_F
    
    dnvec=[dnvec,diffn]; 
    
    Mmaxvec1=[Mmaxvec1,Mmax1];
    
    Fvec1=[Fvec1, Frange];
    Fvec2=[Fvec2, Frange]; 
    Fvec3=[Fvec3, Frange];
    Fvec4=[Fvec4, Frange];
    
    Md= (kB*T)/(2*pi)*log(omg1/(4*pi*fc1)) + c1*sqrt(Frange);
    Mvec1=[Mvec1,Md];
    
    Md= (kB*T)/(2*pi)*log(omg2/(4*pi*fc2)) + c2*sqrt(Frange);
    Mvec2=[Mvec2,Md];
    
    Md= (kB*T)/(2*pi)*log(omg3/(4*pi*fc3)) + c3*sqrt(Frange);
    Mvec3=[Mvec3,Md];
    
    Md= (kB*T)/(2*pi)*log(omg4/(4*pi*fc4)) + c4*sqrt(Frange);
    Mvec4=[Mvec4,Md];
    
    diffM=Mmax1-M3_p;
    diffMvec = [diffMvec,diffM];
    slp = diffM/diffn;
    slpvec=[slpvec;slp];
      
end    



figure
%plot(Fvec1(2:length(Fvec1)) , Mvec1(2:length(Mvec1)), 'b','LineWidth',1,'DisplayName','omg/fc=2'); %2

%hold on
%plot(Fvec2(2:length(Fvec2)) , Mvec2(2:length(Mvec2)), 'r','LineWidth',1,'DisplayName','omg/fc=1'); %1
%hold on
%plot(Fvec3(2:length(Fvec3)) , Mvec3(2:length(Mvec3)), 'g','LineWidth',1,'DisplayName','omg/fc=1/5'); % 1/5
%hold on 
%plot(Fvec4(2:length(Fvec4)) , Mvec4(2:length(Mvec4)), 'y','LineWidth',1,'DisplayName','omg/fc=1/10'); % 1/10
%hold on 
plot(Fvec1(2:length(Fvec1)) , Mpvec1(2:length(Mpvec1)), '--','LineWidth',1,'DisplayName','M_{3_p}'); %2
hold on
plot(Fvec1(2:length(Fvec1)) , Mmaxvec1(2:length(Mmaxvec1)), '.','LineWidth',1,'DisplayName','M_{max}'); %2
xlabel('F (pN)');
ylabel('M_{max} (pNnm)');
legend('show');
title('moment vs force F');

figure
plot(Fvec1(2:length(Fvec1)) , dnvec(2:length(dnvec)), '-','LineWidth',1,'DisplayName','dn');
xlabel('F (pN)');
ylabel('\delta n (turn)');
title('\delta n vs force F');

figure
plot(Fvec1(2:length(Fvec1)),diffMvec(2:length(diffMvec)),'LineWidth',1);
xlabel('F (pN)');
ylabel('M_{max} - M_{3_p}');
title('M_{max} - M_{3_p} vs force F')

figure
plot(Fvec1(2:length(Fvec1)),slpvec(2:length(slpvec)),'LineWidth',1);
xlabel('F (pN)');
ylabel('slope');
title('slope(of stair-step) vs force F')

