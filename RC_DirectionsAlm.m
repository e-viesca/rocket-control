%{
GNU General Public License v2.0
    Copyright (C) 2015  Eugenio Viesca Revuelta, Jos� Carlos Ure�a

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
%                Rocket deviation as a function of angle.
%==========================================================================
% Here we are going to use a loop to plot the deviations caused by the
% Coriolis and centripetal accelerations a function of the different
% launch directions.
% WARNING: Long running time (~15 mins in a 2013 MacBook Air).

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

% air friction
k = 0.001;

% angular rotation speed of the earth
w = [0; 0; 7.27e-5];

% memory allocation for deviation data
dCe = zeros(1,361);
dCo = zeros(1,361);

% we calculate the deviation for each one of the 360�
for I = 0:360
    
    %========
    % ROCKET
    %========
    
    % Initial conditions
    % Almeria: 36� 50' 17'' N   2� 27' 35'' W
    Lat = (90-36.84) * (pi/180);
    Lon = (  -02.46) * (pi/180);
    
    r0 = R * [sin(Lat) * cos(Lon); sin(Lat) * sin(Lon); cos(Lat)];
    
    AlphaV = (90-  30  ) * (pi/180);
    AlphaH = (      I  ) * (pi/180);
    
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
    
    hCe = @(r,v,t) (ag(r,v,t) + af(r,v,t) + aC(r,v,t));
    [rCe,vCe,tCe] = RC_RK2(g,hCe,r0,v0,t0,tN,N);
    
    hCo = @(r,v,t) (ag(r,v,t) + af(r,v,t) + ac(r,v,t));
    [rCo,vCo,tCo] = RC_RK2(g,hCo,r0,v0,t0,tN,N);
        
    %==============
    % CALCULATIONS
    %==============
    
    % we calculate when the rocket lands, then we calculate the actual
    % distance between two shots.
    i=2;
    while norm(r(:,i))>R
        i=i+1;
    end
    
    iCe=2;
    while norm(rCe(:,iCe))>R
        iCe=iCe+1;
    end
    dCe(I+1) = norm((rCe(:,iCe))-(r(:,i)));
    
    iCo=2;
    while norm(rCo(:,iCo))>R
        iCo=iCo+1;
    end
    dCo(I+1) = norm((rCo(:,iCo))-(r(:,i)));
    
end

% we save the data to save us calculate them again in the other two scripts
save('RC_dCedCoAlm.mat','dCe','dCo');

aH = linspace(0,360,361);

figure(1), plot(aH, dCe,'r')
    title('Lanzamiento desde Almer�a. Desviaci�n sin aceleraci�n centr�peta')
    xlabel('Direcci�n de lanzamiento [grados]')
    ylabel('Distancia [m]')
figure(2), plot(aH, dCo,'b')
    title('Lanzamiento desde Almer�a. Desviaci�n sin aceleraci�n de Coriolis')
    xlabel('Direcci�n de lanzamiento [grados]')
    ylabel('Distancia [m]')
