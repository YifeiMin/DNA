%% set myfun_transition1=1 to solve for Mcric and L_p1, thus get critical n


global L Kb Kt T b kB F bp l_bp lB cr Wr dn a kd v eta nlp;



format long

x0_s=[24.6,L*0.08]; %initial searching value
[x_s,fval_s] = fsolve(@myfun_transition1,x0_s);  % Call solver to solve for Mcric and Lp

Mcric_s= x_s(1);
Lp0_s = x_s(2);
K_s= ((Kb*F-((Mcric_s)^2)/4)^(1/2))/(kB*T);  % K in eqn(6)
ncric_s = ((Mcric_s*L)/(2*pi))*((1/Kt)+(1/(4*Kb*K_s)));





