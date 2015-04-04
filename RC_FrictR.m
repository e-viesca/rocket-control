%{
GNU General Public License v2.0
    Copyright (C) 2015  Eugenio Viesca Revuelta, José Carlos Ureña

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
%}

%==========================================================================
%                             ROCKET CONTROL
%             Travel distance as a function of air friction.
%==========================================================================
% Here we are going to use a loop to plot the travel distance as a function
% of a range of possible coefficients of air friction.
% ATTENTION: Prolonged running time (~2 mins).

tic;

clear all
close all
clc

%=======
% EARTH
%=======

% Gravity
G = 6.667e-11;
R = 6371009;
M = 5.9736e24;

% angular rotation speed of the earth
w = [0; 0; 7.27e-5];

% loop iterations
n = 100;

% air friction
k0 = 0.0006;
k  = k0;
kN = 0.05;
dk = (kN-k0)/n;

% memory allocation for distance data
dR = zeros(1,n);

% we calculate the distance for each one of the n-different k's
for I = 1:n;
    
    %========
    % ROCKET
    %========
    
    % Initial conditions
    % Almeria: 36º 50' 17'' N   2º 27' 35'' W
    Lat = (90-36.84) * (pi/180);
    Lon = (  -02.46) * (pi/180);
    
    r0 = R * [sin(Lat) * cos(Lon); sin(Lat) * sin(Lon); cos(Lat)];
    
    AlphaV = (90-  45  ) * (pi/180);
    AlphaH = (     45  ) * (pi/180);
    
    e1 = r0 ./ norm(r0);
    e2 = R * [ cos(Lat) * cos(Lon); cos(Lat) * sin(Lon); -sin(Lat)];
    e2 = e2./norm(e2);
    e3 = R * [-sin(Lat) * sin(Lon); sin(Lat) * cos(Lon);     0    ];
    e3 = e3./norm(e3);
    
    v0 = 7000 * [sin(AlphaV) * cos(AlphaH); sin(AlphaV) * sin(AlphaH); cos(AlphaV)];
    v0 = [e3, e2, e1] * v0;
    
    %=============
    % INTEGRATION
    %=============
    t0 = 0;
    tN = 1.1*(2*norm(v0)/(G*M/R^2));
    N = 10000;
    
    %===============
    % ACCELERATIONS
    %===============
    % gravity
    ag = @(r,v,t)  ( -((G*M)./norm(r)^3 .* r) );
    % air friction
    af = @(r,v,t)  ( -k .* v);
    % centripetal force
    ac = @(r,v,t)  ( cross(w, cross(w, r)) );
    % Coriolis force
    aC =  @(r,v,t) (-cross(2.*w, v) );
    
    
    %========================
    % DIFFERENTIAL EQUATIONS
    %========================
    g = @(r,v,t) (v);
    h = @(r,v,t) (ag(r,v,t) + af(r,v,t) + ac(r,v,t) + aC(r,v,t));
    
    [r,v,t] = RC_RK2(g,h,r0,v0,t0,tN,N);
    
    
    %==============
    % CALCULATIONS
    %==============
    
    % we calculate when does the rocket land, then we calculate the
    % distance to the launching site.
    i=2;
    while norm(r(:,i))>R
        i=i+1;
    end
    
    dR(I) = norm((r(:,1))-(r(:,i)));
    
    k = k + dk;
    
end

K = linspace(k0,k,n);

figure(1), plot(K, dR,'r')
    title('Distancia recorrida en función del coeficiente de rozamiento')
    xlabel('k')
    ylabel('Distancia [m]')

toc