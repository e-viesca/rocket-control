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

close all
clear all
clc

%=======
% EARTH
%=======

% Gravity
G = 6.667e-11;
R = 6371000;
M = 5.9736e24;

% air friction
k = 0.001; % to be played with

% angular rotation speed of the earth
w = [0; 0; 7.27e-5];

%========
% ROCKET
%========

% initial conditions
% Almeria: N36º 50' 17'' W2º 27' 35''
Lat = (90-36.84) * (pi/180); 
Lon = (2.46) * (pi/180);
r0 = [R * sin(Lat) * cos(Lon), R * sin(Lat) * sin(Lon), R * cos(Lat)];

AlphaV = (90-45) * (pi/180);
AlphaH = (-45) * (pi/180);

e3 = [sin(AlphaV) * cos(AlphaH)];
e2 = [sin(AlphaV) * sin(AlphaH)];
e1 = [cos(AlphaV)];

v0 = 5000 * [e1 , e2, e3];
% v0 = (...);

%=============
% INTEGRATION
%=============
t0 = 0;
tN = 2*norm(v0)/(G*M/R^2);
N = 1000;

%===============
% ACCELERATIONS
%===============
% gravity
ag = @(r,v,t) ( -((G*M)./(r.^2)) );
% air friction
af = @(r,v,t) ( -k * v);
% centripetal force
ac = @(r,v,t) ( cross(w, cross(w, r)) );
% Coriolis force
a =  @(r,v,t) ( cross(2.*w, v) );

%{
% propulsion (NO NO for now)
ap = @(r,v,t) ();
dmdt = @(t) ();
%}

%========================
% DIFFERENTIAL EQUATIONS
%========================
g = @(r,v,t) ();
h = @(r,v,t) ();

[r,v,t] = rocketRK2(g,h,r0,v0,t0,tN,N);

%===============
% VISUALISATION
%===============
figure(1)
hold on

% plot starting point
plot3(r0(1),r0(2),r0(3),'rx','MarkerSize',20,'LineWidth',4)

% plot trajectory
%plot3(r(1,:),r(2,:),r(3,:),'.','LineWidth',20)

% plot Earth
[x,y,z] = sphere(30);
mesh(R*x,R*y,R*z)
colormap([0 0 0])
view(-180,20)
axis equal