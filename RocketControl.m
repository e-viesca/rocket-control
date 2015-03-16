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

close all
clear all
clc

%=======
% EARTH
%=======

% Gravity
G = 6.667e-11;
R = 6371009;
M = 5.9736e24;

% air friction
k = 0.001; % to be played with

% angular rotation speed of the earth
w = [0; 0; 0]; %7.27e-5];

%========
% ROCKET
%========

% Initial conditions
% Almeria: N36� 50' 17'' W2� 27' 35''
Lat = (90-36.84) * (pi/180); 
Lon = (2.46) * (pi/180);
%Lat = (0) * (pi/180); 
%Lon = (0) * (pi/180);
r0 = [R * sin(Lat) * cos(Lon); R * sin(Lat) * sin(Lon); R * cos(Lat)];

AlphaV = (90-  45  ) * (pi/180);
AlphaH = (     45  ) * (pi/180);

e1 = r0 ./ norm(r0);
e2 = R * [ cos(Lat) * cos(Lon); cos(Lat) * sin(Lon); -sin(Lat)];
e2 = e2./norm(e2);
e3 = R * [-sin(Lat) * sin(Lon); sin(Lat) * cos(Lon);     0    ];
e3 = e3./norm(e3);

v0 = 1 * [sin(AlphaV) * cos(AlphaH); sin(AlphaV) * sin(AlphaH); cos(AlphaV)];
%%%%%v0 = [e1, e2, e3] * v0;

%=============
% INTEGRATION
%=============
t0 = 0;
tN = 2*norm(v0)/(G*M/R^2);
N = 10000;

%===============
% ACCELERATIONS
%===============
% gravity
ag = @(r,v,t) ( -((G*M)./norm(r)^3 .* r) );
% air friction
af = @(r,v,t) ( -k .* v);
% centripetal force
ac = @(r,v,t) ( cross(w, cross(w, r)) );
% Coriolis force
aC =  @(r,v,t) ( cross(2.*w, v) );

%{
% propulsion (NO NO for now)
ap = @(r,v,t) ();
dmdt = @(t) ();
%}

%========================
% DIFFERENTIAL EQUATIONS
%========================
g = @(r,v,t) (v);
h = @(r,v,t) (ag(r,v,t) + af(r,v,t) + ac(r,v,t) + aC(r,v,t));

[r,v,t] = rocketRK2(g,h,r0,v0,t0,tN,N);
%===============
% VISUALISATION
%===============
quiver3(0,0,0,v0(1),v0(2),v0(3))
%%
figure(1)
hold on

% plot starting point
plot3(r0(1),r0(2),r0(3),'rx','MarkerSize',20,'LineWidth',3)

% plot trajectory
plot3(r(1,:),r(2,:),r(3,:),'.','LineWidth',20)

% plot Earth
[x,y,z] = sphere(30);
mesh(R*x,R*y,R*z)
colormap([0 0 0])
view(-180,20)
axis equal