%%algebraic equations including undulation dr  
function matrix = myfun_ubbink(x)
global Kb T kB F lB cr v kd eta
 % Marko and Siggia Electrostatics + Entropy 150 mm
u = x(1); %theta
r = x(3); 
M3 = x(2); 
dr = x(4); 

%% cambiar esto
A = (Kb/kB/T)^(1/3); 
entropy = cr*kB*T/(A*(dr^(2/3))) + eta*kB*T/(A*((r*cot(u)))^(2/3)); %Uconf %eqn(19)
dentropydr =  eta*(-2/3)* kB*T/(A*(r*cot(u))^(5/3))*(cot(u)); %dr
deddr = cr*(-2/3)* kB*T/(A*dr^(5/3)); % d dr

dentropydtheta = eta*(2/3)*(kB*T*r)*(1+(cot(u))^2)/(A*(r*cot(u))^(5/3)); 

%%Poisson Boltzman

C = (1/2)*kB*T*v^2*lB*(pi/kd/r)^(1/2);
g = 1 + 0.83*(tan(u))^2+0.86*(tan(u))^4;

Upb = C*g*exp(2*kd^2*dr^2-2*kd*r) + entropy ;
dupb_dr= C*g*exp(2*kd^2*dr^2-2*kd*r)*4*kd^2*dr + deddr; % d dr
DUpb_r = -2*C*kd*exp(2*kd^2*dr^2-2*kd*r)*g + dentropydr; %d r
DUpb_t = C*exp(2*kd^2*dr^2-2*kd*r)*(1.66*tan(u)*(1+(tan(u))^2) +3.44*(tan(u))^3*(1+(tan(u))^2) ) + dentropydtheta;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% theta_o_dot


%%Additions
K = ((Kb*F-(x(2)^2)/4)^(1/2))/(kB*T);
Gflu = ((kB*T)^2*K/Kb)*(1-1/(4*K)-1/(64*K^2));


%% Lkhat = x(2)*(x(2)*k*T/(4*Kb*(Kb*F)^(1/2)));
Lkhat = x(2)*(x(2)*kB*T/(4*Kb*(Kb*F-(x(2)^2)/4)^(1/2))); %M^2/4KbK

%%Eq2 = Kb*(-2*sin(theta)^3*cos(theta))/(r^2)+M3*cos(2*theta)/r;
%%Eq1= Kb/2*(sin(theta)^4/r^2 ) + Upb - M3*sin(2*theta)/2/r + F

% Explanation: Lagrange multiplier is replace by \lambda=2pi*M3, which is given by setting dU/dM3 = 0 
% Explanation: Lp doesn't appear, since LK num will later be represented as a fct of Lp

Eq1= Kb/2*(sin(x(1))^4/r^2 ) + Upb - x(2)*sin(2*x(1))/2/r + F - Gflu + 0*Lkhat ; % d Lp
Eq2 = Kb*(-2*sin(x(1))^3*cos(x(1)))/(r^2) + x(2)*cos(2*x(1))/r - DUpb_t; %d theta
Eq3 = -Kb*sin(x(1))^4/r^3 + DUpb_r + x(2)*sin(2*x(1))/(2*r^2); %d r
Eq4 = dupb_dr;  % d dr
matrix = [Eq1;Eq2;Eq3;Eq4];

end