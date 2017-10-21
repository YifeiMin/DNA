%% energy of the straight state
global L Kb Kt T kB F ;
 
%% energy _ straight regime with fluctuation
n= Lk; % linking number 
fener=@(x) (x/(2*(pi)))*L*((1/Kt)+(1/(4*((sqrt(Kb*F - ((x^2)/4)))/(kB*T))*Kb)))-n;
M_s = fzero(fener,0); %pNnm %moment of the straight regime as a function of linking number

K_s = (sqrt(Kb*F - ((M_s^2)/4)))/(kB*T);
Gflu = (((kB*T)^2)/Kb)*K_s*(1-(1/(4*K_s))-(1/(64*(K_s^2))));
V_straight = (-F + ((M_s^2)/(2*Kt))+Gflu+((M_s^2)/(4*K_s*Kb)) )*L;
