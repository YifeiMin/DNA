global  T kB Kb  F L Kt;

M0 = 25; % initial search value pNnm
[x,fval] = fsolve(@findMcrit1, M0);