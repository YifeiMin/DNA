function dummyeq= findMcrit1(x)
global Kb Kt L kB T F 

dummyeq= 22.36-(x/(2*pi))*L*(1/Kt + 1/(4*Kb*(sqrt(Kb*F - ((x^2)/4))/(kB*T))));
end

