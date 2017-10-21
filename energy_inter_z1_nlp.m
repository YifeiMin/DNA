% This script computes the linking number of the rod
% z_2 and z_3 are set to be 1
%z_1 can be adjusted
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remark: all the values in this script are scaled (including curvature twisting etc.), unless specified 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% R phi Z are centerline cylindrical coordinates
% s is the arclength with scaling
%set value of z1 first
global L Kb Kt T b kB F z1 l_bp nlp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%set unit for length, and torque


lambda=sqrt(Kb/F); % nm  length unit   
Torqueunit = sqrt(Kb*F);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% set value for tau
tau = sqrt((1+z1)/(1-z1));
%%

%Euler angle as a function of arclength, actual form
syms s;

varphi(s)=atan((1/(tau))*((sinh((s/lambda)))/(cosh((s/lambda))))) + tau*(s/lambda);
z(s)=1-(2/(1+tau^2))*(sech((s/lambda)))^2;
psi(s)=atan((1/tau)*((sinh((s/lambda)))/(cosh((s/lambda))))) + (3 - (2/b))*tau*(s/lambda);

theta(s) = acos(z);
%%

% twist vertor as a function of arclength, actual form
kappa1(s)=diff(theta(s))*sin(psi(s))-diff(varphi(s))*sin(theta(s))*cos(psi(s));
kappa2(s)=diff(theta(s))*cos(psi(s))+diff(varphi(s))*sin(theta(s))*sin(psi);
kappa3(s)=diff(varphi(s))*z(s)+diff(psi(s));

kap1=@(s) eval(kappa1);
kap2=@(s) eval(kappa2);
kap3=@(s) eval(kappa3);

t = ((-L/2):0.01:(L/2)); %step size

k1=kap1(t);
k2=kap2(t);
k3=kap3(t);

%plot twist vector
% figure;
% plot(k1,t,'color','black','linewidth',1); 
% hold on;
% plot(k2,t,'color','b','linewidth',1); 
% hold on;
% plot(k3,t,'color','g','linewidth',1); 
%%

fR =lambda*(2/(1+(tau)^2))*sech(s/lambda) ;

fphi = tau*(s/lambda) - pi/2;

%transform into cartisian coordinates, actual form
fX = fR*cos(fphi);
fY = fR*sin(fphi);
fZ = s- lambda*(2 / (1+tau^2))*(sinh(s/lambda)/cosh(s/lambda)) ;


Xs = @(s) eval(fX);
Ys = @(s) eval(fY);
Zs = @(s) eval(fZ);

%plot centerline 

X=Xs(t);
Y=Ys(t);
Z=Zs(t);

% figure;
% plot(t,Zs(t));

% figure;
% plot3(X,Y,Z, 'linewidth',2);
% axis([-(L/2) (L/2) -(L/2) (L/2) -(L/2) (L/2) ])

M3actual = sqrt((1+z1)*(1+1)*(Kb*F)); %pNnm

%calculate linking number 
Tw=(0.5/(pi))*integral(kap3,-(L/2),(L/2)); %Twist

Wrf= -((z(s)-1)*diff(varphi(s)));
Wrf1=@(s) eval(Wrf);
Wr=(0.5/(pi))*integral(Wrf1,-(L/2),(L/2)) ; %writhe

Lk=Tw+Wr*nlp; %linking number 
Ztotal=Zs(L/2)-Zs(-L/2);  % vertical displacement of twisted DNA

difZ = 1-(2/(1+tau^2))/((cosh(s/lambda))^2);
tolo=0.03; % tolerance
difZ1= @(s) eval(difZ)-(1-tolo); %set to zero to solve for lo

Lo = 2*abs(fzero( @(s)difZ1(s),(-L/(l_bp/500)))); %nm

Lt = L - nlp* Lo; %nm

K= sqrt(Kb*F - ((M3actual^2)/4))/(kB*T); 

Etwist_t = (M3actual^2)*Lt/(2*Kt); %pNnm
Ewrith_t = (M3actual^2)*Lt/(4*Kb*K); %pNnm
Gflu= (((kB*T)^2)/Kb)*K*(1-(1/(4*K))-(1/(64*(K^2))));
Eflu_t = Gflu*Lt;
Vtail=Etwist_t + Ewrith_t + Eflu_t ; %w/o potential of tension

Velas_o=nlp*(F*Lo+((M3actual^2)*Lo)/(2*Kt));

V_3=Vtail+Velas_o - F*((Zs(L/2)-Zs(-L/2))-(nlp-1)*( Lo-(Zs(Lo/2)-Zs(-Lo/2))  )); %pNnm  energy of the 3rd state


