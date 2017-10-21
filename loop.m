
function matrix = loop(x,entry)
global  T kB Kb Kt L Wr dn a v kd lB cr eta


I = entry(1);
M3 = entry(2);
u = entry(3);
r = entry(4);
dr = entry(5);


A = (Kb/kB/T)^(1/3);
entropy = cr*kB*T/(A*(dr^(2/3))) + eta*kB*T/(A*((r*cot(u)))^(2/3));
C = (1/2)*kB*T*v^2*lB*(pi/kd/r)^(1/2);
g = 1 + 0.83*(tan(u))^2+0.86*(tan(u))^4;

Upb = C*g*exp(2*kd^2*dr^2-2*kd*r) + entropy ; 

%% Marko and Siggia Electrostatics
dM= x(1); 
Lp = x(2); 

M = M3 +dM; 

%%Additions
P3 = ((Kb.*I-(M3.^2)/4)^(1/2))./(kB.*T);  % K in eqn(6)

Gflu3 = ((kB.*T)^2.*P3./Kb).*(1-1./(4.*P3)-1./(64.*P3.^2));

Lkhat3 = M3.^2./(4.*Kb.*P3); 

%%Additions
P = ((Kb.*I-(M.^2)/4)^(1/2))./(kB.*T); % K in eqn(6)

Gflu = ((kB.*T)^2.*P./Kb).*(1-1./(4.*P)-1./(64.*P^2));

Lkhat = M.^2./(4.*Kb.*P);

Lo = 4*(Kb/I)^(1/2);

kappa = ((sin(u))^2)/r;

e1 = L*(- I + (M)^2/(2*Kt) + Gflu + Lkhat ); %energy of straight regime
e2 = L*(- I + (M3)^2/(2*Kt) + Gflu3 + Lkhat3 ) + Lo*( 2*I + a*( - Gflu3 - Lkhat3) ) + (Lp)*(I + Upb - Gflu3 - Lkhat3 )+  1/2*(kappa^2*Kb*Lp) ; %energy of plectonemic regime
eq1 = e1-e2;

e3 = M*L/(2*pi)*(1/Kt + 1/(4*Kb*P)); % Linking num of straight
e4 = M3*L/(2*pi)*(1/Kt + 1/(4*Kb*P3)) - (a*Lo+Lp)*M3/(2*pi)*(1/(4*Kb*P3)) + (Wr-dn) + sin(2*u)/(4*pi*r)*Lp; %link of plectoneme

eq2 = e3-e4 ; 

matrix = [eq1;eq2];

end
