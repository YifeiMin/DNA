%% this script find transition point(i.e. critical link and moment)

function matrix = myfun_transition1(x)

Mcric = x(1);
Lp0=x(2);

global L Kb Kt T b kB F bp l_bp  lB cr Wr dn a kd v eta ;

run_ubbink_1  %%   (Amoment)co = 150 mM  Marko and Siggia Electrostatics and Entropy %% solve for configuration of plectoneme

%% parameter of plectoneme regime
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

%% parameter of straight regime
K_s= ((Kb*F-((Mcric)^2)/4)^(1/2))/(kB*T);  % K in eqn(6)
Gflu_s= ((kB*T)^2*(K_s)/Kb)*(1-1/(4*(K_s))-1/(64*(K_s)^2));


%% set following to be equal
V_p1 = L*(- F + (M3_p)^2/(2*Kt) + Gflu_p + Lkhat_p ) + Lo_p*( 2*F + a*( - Gflu_p - Lkhat_p) ) + (Lp0)*(F + Upb_p - Gflu_p - Lkhat_p )+  1/2*(kappa_p^2*Kb*Lp0) ; %energy of plectonemic regime as fct of Lp1

n_p1 = M3_p*L/(2*pi)*(1/Kt + 1/(4*Kb*K_p)) - (a*Lo_p+Lp0)*M3_p/(2*pi)*(1/(4*Kb*K_p)) + (Wr-dn) + sin(2*theta_p)/(4*pi*r_p)*Lp0; %link of plectoneme as fct of Lp1

V_s1=(-F+(Mcric^2)/(2*Kt)+Gflu_s+(Mcric^2)/(4*Kb*K_s))*L;

n_s1=((Mcric*L)/(2*pi))*((1/Kt)+(1/(4*Kb*K_s)));

eq1=V_p1-V_s1;
eq2=n_p1-n_s1;

matrix=[eq1;eq2];

end