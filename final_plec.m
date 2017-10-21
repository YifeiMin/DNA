clc;
clear all;
global  T kB Kb Kt lB cr L Wr dn a kd v eta


format long
T = 296.5; % Kelvin 296.5 for forth, 300 for them
kB = 1.3806503*10^(-2); %pNnm K^-1
Kb = 50*kB*T; %pNnm^2 45 for 320 46 for for
Kt = 95*kB*T;%pNnm^2 100 for 320
lB= 0.715; %nm
bp = 0.34; %nm %length of a base pair
l_bp = 2200;
L = l_bp*bp;
cr = 1/(2^(8/3)); % parameter control entropy
eta = cr/(pi)^(2/3);
Wr = 1; % parameter control size loop
dn = 0; % delta n at the transiton 0.7 for short and 1 for long 
a=1;  % parameter control fluctuations
w=0; % parameter controls moment assumption in the loop. w =0 means no moment in
% the loop (Kulic)
y= 1;
v = 0.72*8.231260437 ;  % 0.72 47.7 ; %8.06 ;%% 1.97; % nm^-1
kd =1/0.784874795;%1/0.351; %1/0.8; %%1/3.07; %mM %reciprocal of lambda_D, eqn(18)


%%2.2 kbp variables
run_hamiltonian_ubbink %%   (Amoment)co = 150 mM  Marko and Siggia Electrostatics and Entropy 

F =[Fext1,Fext2]; 

%% Critical torque from loop analysis 
Fs = F; %%[0.2:0.1:4]; % vector of different tension forces
%% Definitions of Moments according to salt concentration
Mvector=[moment1;moment2]; % Moments in the plectoneme for 320 mM
thetavector = [theta1;theta2];
rvector = [radius1;radius2];
rrvector = [rr1;rr2];


%Kulic + Coyne expressionsm + thermal fluctuations
x=[0,0];
dM = [];
dLp =[];

for i=1:length(F) 
format long
M1 = x; 
entry = [F(i),Mvector(i),thetavector(i),rvector(i),rrvector(i)];
[x,fval] = fsolve(@(x)loop(x,entry),M1);  % find deltaM and Lp
dM = [dM;x(1)];
dLp = [dLp;x(2)];
end



figure
plot(F,dM, 'b-','LineWidth',2)
hold on
xlabel('F[pN]','fontsize',24,'fontweight','b'); ylabel('M_3[pNnm]','fontsize',24,'fontweight','b');
axis tight

lo = y*4*(Kb./Fs').^(1/2); % loop length
Mcrik2 = dM+Mvector;

%%
%% PLOTS 1.9 AND 10.9 kbp
Fexp = [1:0.5:3.5]-0.025;
Mexp =[17.25 20.5 22.77 25.66 27.25 28.48];


figure
plot([Fext1,Fext2],Mvector, 'b-','LineWidth',2) % Moments in the plectoneme for 150 mM
hold on
plot(Fexp,Mexp, 'blacko','LineWidth',5)
hold on
plot(Fs,Mcrik2, 'blue','LineWidth',2) % 2.2 kbp
xlabel('F[pN]','fontsize',24,'fontweight','b'); ylabel('M_3[pNnm]','fontsize',24,'fontweight','b');
axis tight



%% 320 mM analysis
l_bp = 2200;
Length1 = l_bp*bp;
link1 = Mcrik2.* (Length1*(1/Kt + kB*T./(4*Kb*(Kb*Fs' - Mcrik2.^(2)/4).^0.5))/(2*pi)); % link complete independent of salt concentration
link2 = Mvector.* (Length1*(1/Kt + kB*T./(4*Kb*(Kb*F' - Mvector.^(2)/4).^0.5))/(2*pi)); % link for 320 mM plectoneme


n_star = [5.25 8.08 9.7];
F_star =[1 2 3];

figure
plot(F,link1, 'r-','LineWidth',2) % McriK2 complete
hold on
plot(F,link2, 'b-','LineWidth',2) % plectoneme
hold on
plot(F_star,n_star, 'ro','LineWidth',10) % experimental
title('Critical Link number vs F', 'fontsize',16,'fontweight','b');
xlabel('F[pN]','fontsize',14,'fontweight','b'); ylabel('n*','fontsize',14,'fontweight','b');
axis tight



dn2 = dn;

r = rvector;
theta = thetavector;
P1 = ((Kb.*F'-(Mvector.^2)/4).^(1/2))./(kB.*T);

lp=(link1-link2-(Wr-dn) + a*Mvector.*lo./(8*pi*Kb*P1)).*(sin(2*theta)./(4*pi*r)- Mvector./(8*pi*Kb*P1)).^(-1);

epsilon1 = ( Kb.*F'/(kB.*T)^2 - Mvector.^2./(2.*kB.*T)^2 - 1/32).^0.5;
row1= 1-1/2.*(1./epsilon1) + Kb.*kB.*T./(Length1.*(Kb.*F'-Mvector.^2/4)) ; % plectoneme

epsilon2 = ( Kb.*F'/(kB.*T)^2 - Mcrik2.^2./(2.*kB.*T)^2 - 1/32).^0.5;
row2= 1-1/2.*(1./epsilon2) + Kb.*kB.*T./(Length1.*(Kb.*F'-Mcrik2.^2/4)) ; % with torque expression complete nm

jump = row1.*(lp+ lo) + (row2-row1).*(Length1);



%%EXPERIMENTAL DATA 150 mM
j_exp =[110.5; 96.3; 76.81; 74.235; 65.37; 62.31];%2.2 kbp
U_error = [(122.2-110.5) ;(101.93-96.3); (80.4-76.81); (77.78-74.235); (67.3-65.37); (65-62.31)];
D_error = [(113.8-111.6); (98.23-90.34); (88.41-82.45); (79.22-73.91); (71.98-69.57)];
f_e = [0.97; 1.46; 1.94; 2.43; 2.92; 3.4];
j_exp2 =[111.6; 98.23; 82.45; 73.91; 69.57];%4.2 kbp
f_e2 = [0.93; 1.43; 1.92; 2.40; 2.88];


figure
plot(F,jump,'b-','LineWidth',2)
hold on
plot(f_e,j_exp,'bo','LineWidth',8) %2.2 kbp
hold on
errorbar(f_e,j_exp,U_error,U_error,'ro','LineWidth',2) %2.2 kbp
hold on
errorbar(f_e2,j_exp2,D_error,D_error,'bo','LineWidth',2) %2.2 kbp
hold on
title('Jump in the extension vs F', 'fontsize',16,'fontweight','b');
xlabel('F[pN]','fontsize',14,'fontweight','b'); ylabel('n*','fontsize',14,'fontweight','b');
axis tight

figure
plot(F,link1-link2 + dn,'b-','LineWidth',2)
title('deltaN vs F', 'fontsize',16,'fontweight','b');
xlabel('F[pN]','fontsize',14,'fontweight','b'); ylabel('n*','fontsize',14,'fontweight','b');
axis tight


P1 = ((Kb.*F'-(Mvector.^2)/4).^(1/2))./(kB.*T);
Ghat1 = ((kB.*T)^2.*P1./Kb).*(1-1./(4.*P1)-1./(64.*P1.^2));
Lkhat1 = Mvector.^2./(4.*Kb.*P1);
kappa = ((sin(thetavector)).^2)./rvector;

A = (Kb/kB/T)^(1/3);
rrvector=[rr1;rr2];
entropy = cr*kB*T./(A*(rrvector.^(2/3))) + eta*kB*T./(A*((rvector.*cot(thetavector))).^(2/3));
C = (1/2)*kB*T*v^2*lB*(pi/kd./rvector).^(1/2);
g = 1 + 0.83*(tan(thetavector)).^2+0.86*(tan(thetavector)).^4;

Upb = C.*g.*exp(2*kd^2*rrvector.^2-2*kd*rvector) + entropy ;
slope = [slope1;slope2];
Uo = (2*F' - Ghat1 - Lkhat1).*lo;
U1 = (F' - Ghat1 - Lkhat1 + 1/2*(kappa.^2*Kb) + Upb).*slope./row1;

%%
%%
%%
%%
%%


figure
plot(F,Uo,'b-','LineWidth',2)
hold on
plot(F,U1,'r-','LineWidth',2)
hold on
plot(F,U1+kB*T,'magenta-','LineWidth',0.5)
hold on
plot(F,U1-kB*T,'magenta-','LineWidth',0.5)
title('Energy per unit lenght', 'fontsize',16,'fontweight','b');
xlabel('F[pN]','fontsize',14,'fontweight','b'); ylabel('n*','fontsize',14,'fontweight','b');
axis tight

F_experimental =  [1:0.5:3.5]-0.025;
p_experimental = [49.28 38.97 35.43 32.2 30.7 29.47];
F_experimental2 =  [1:0.5:3]+0.02;
p_experimental2 = [51.49 41.18 38.45 33.61 29.6];


figure
plot(F',slope,'b-','LineWidth',2)
hold on
plot(F_experimental,p_experimental,'ro','LineWidth',8)
hold on
plot(F_experimental2,p_experimental2,'ro','LineWidth',8)
hold on
plot(F',row1.*lo,'b-','LineWidth',2)
axis tight


[Uo-U1+kB*T./lo,Uo-U1-kB*T./lo, Uo-U1,F']
[F',[radius1;radius2]]
