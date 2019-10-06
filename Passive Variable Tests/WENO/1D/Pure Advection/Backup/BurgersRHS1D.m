function [du] = BurgersRHS1D(x,u,h,k)

% function [du] = BurgersRHS1D(x,u,h,k)
% Purpose  : Evaluate the RHS of Burgers equations using 
% an ENO/WENO scheme

N = length(x);
du = zeros(N,1);

xl = zeros(N+2*(k-1),1);
ul = zeros(N+2*(k-1),1);

% Extend x and u
if (k==1)
    xl = x; ul = u;
else
  for q=1:(k-1)
    ul(k-q) = u(N-(q-1)); xl(k-q) = -q*h;
    ul(N+k+q-1) = u(q); xl(N+k+q-1) = 1 + (q-1)*h;
  end;
  ul(k:(N+k-1)) = u(1:N);
  xl(k:(N+k-1)) = x(1:N);
end

% define internal left and right interface values
uml = zeros(N,1); umr = zeros(N,1);
uloc = zeros(2*k-1,1); xloc = zeros(2*k-1,1);

for i=1:N
  xloc = xl(i:(i+2*(k-1))); uloc = ul(i:(i+2*(k-1)));  
  % [uml(i),umr(i)] = ENO(xloc,uloc,k);
   [uml(i),umr(i)] = WENO(xloc,uloc,k);
end;

% Compute residual at x=0
du(1) = -( numflux(umr(1),uml(2)) - numflux(umr(N),uml(1)))/h;

% Compute residual at x=1-h;
du(N) = -( numflux(umr(N),uml(1)) - numflux(umr(N-1),uml(N)))/h;

% Compute residual every where else
for i=2:(N-1);
   du(i) = -( numflux(umr(i),uml(i+1)) - numflux(umr(i-1),uml(i)))/h;
end

return