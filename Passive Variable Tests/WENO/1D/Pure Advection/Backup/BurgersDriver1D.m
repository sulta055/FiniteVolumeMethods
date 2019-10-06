% Driver script for solving the 1D Burgers equations 
% with developing shock. Periodic boundary conditions

% Solver is ENO or WENO method of order k-1/2k-1 in space and a 3rd order
% SSP RK in time

% N: Number of grid points
% k: width of stencil (k=1,2,3)

N = 100;
k = 2;
CFL=0.1;

%
h = 1/N;
x = [0:h:1-h]';

% Set initial conditions as cell average
u = cellave(x,0,h);

% Solve Problem
FinalTime = 1.0;
[u,time] = Burgers1D(x,u,N,k,CFL,FinalTime,1);

plot(x,u,x,cellave(x,FinalTime,h))

h*sum(u-cellave(x,FinalTime,h))

h*sum(abs(u-cellave(x,FinalTime,h)))

