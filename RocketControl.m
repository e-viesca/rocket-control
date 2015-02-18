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
%%%%% r0 = (...);

AlphaV = (90-45) * (pi/180);
AlphaH = (-45) * (pi/180);
%{
e3 = (...)
e2 = (...)
e1 = (...)

v0 = 5000 * (...)
v0 = (...)
%}

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
ag = @(r,v,t) ();
% air friction
af = @(r,v,t) ();
% centripetal force
ac = @(r,v,t) ();
% Coriolis force
a = @(r,v,t) ();

% propulsion (NO NO by the moment)
ap = @(r,v,t) ();
dmdt = @(t) ();

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
plot3(r0(1),r0(2),r0(3),'gx','MarkerSize',15)

% plot trajectory
plot3(r(1,:),r(2,:),r(3,:),'.','LineWidth',20)

% plot Earth
[x,y,z] = sphere(30);
mesh(R*x,R*y,R*z)
colormap([0 0 0])
view(-180,20)
axis equal