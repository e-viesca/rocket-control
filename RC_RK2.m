function [ r,v,t ] = RC_RK2( g, h, r0, v0, t0, tN, N )
% RC_RK2 Modified 2nd order Runge-Kutta function for two coupled ODE.


% initialisation block
%----------------------
% step size
dt   = (tN-t0)/(N-1);

% arrays
t    = linspace(t0,tN,N);
r    = [zeros(1,N);zeros(1,N);zeros(1,N)];
v    = [zeros(1,N);zeros(1,N);zeros(1,N)];
%t(1) = t0; this is automatically initialized by linspace()
r(:,1) = r0;
v(:,1) = v0;


% integration block
%-------------------
for i=1:N-1
    
    % mid-point values
    tmid = t(i) + dt/2;
    xmid = r(:,i) + g(r(:,i),v(:,i),t(i)) .* dt/2;
    vmid = v(:,i) + h(r(:,i),v(:,i),t(i)) .* dt/2;
    
    % actual step
    r(:,i+1) = r(:,i) + g(xmid,vmid,tmid).*dt;
    v(:,i+1) = v(:,i) + h(xmid,vmid,tmid).*dt;
        
end

