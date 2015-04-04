function [ r,v,t ] = RC_predictionRK2( g, h, r0, rN, t0, tN, N )
% RC_predictionRK2 Modified 2nd order Runge-Kutta function for two coupled
% ODE.
%==========================================================================
%                              DISCONTINUED
%      Failed attempt to build an inverse Runge Kutta for predictor
%==========================================================================


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
r(:,N) = rN;


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

