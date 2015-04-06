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
%               v(t) for three different drag coefficients
%==========================================================================

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

% air friction 1
k1 = 0.0006;
% gravity
ag = @(r1,v1,t)  ( -((G*M)./norm(r1)^3 .* r1) );
% air friction
af = @(r1,v1,t)  ( -k1 .* v1);
% centripetal force
ac = @(r1,v1,t)  ( cross(w, cross(w, r1)) );
% Coriolis force
aC =  @(r1,v1,t) (-cross(2.*w, v1) );
% ODEs
g = @(r1,v1,t) (v1);
h = @(r1,v1,t) (ag(r1,v1,t) + af(r1,v1,t) + ac(r1,v1,t) + aC(r1,v1,t));

[r1,v1,t] = RC_RK2(g,h,r0,v0,t0,tN,N);


% air friction 2
k2 = 0.001;
% gravity
ag = @(r2,v2,t)  ( -((G*M)./norm(r2)^3 .* r2) );
% air friction
af = @(r2,v2,t)  ( -k2 .* v2);
% centripetal force
ac = @(r2,v2,t)  ( cross(w, cross(w, r2)) );
% Coriolis force
aC =  @(r2,v2,t) (-cross(2.*w, v2) );
% ODEs
g = @(r2,v2,t) (v2);
h = @(r2,v2,t) (ag(r2,v2,t) + af(r2,v2,2) + ac(r2,v2,t) + aC(r2,v2,t));

[r2,v2,t] = RC_RK2(g,h,r0,v0,t0,tN,N);


% air friction 1
k3 = 0.002;
% gravity
ag = @(r3,v3,t)  ( -((G*M)./norm(r3)^3 .* r3) );
% air friction
af = @(r3,v3,t)  ( -k3 .* v3);
% centripetal force
ac = @(r3,v3,t)  ( cross(w, cross(w, r3)) );
% Coriolis force
aC =  @(r3,v3,t) (-cross(2.*w, v3) );
% ODEs
g = @(r3,v3,t) (v3);
h = @(r3,v3,t) (ag(r3,v3,t) + af(r3,v3,t) + ac(r3,v3,t) + aC(r3,v3,t));

[r3,v3,t] = RC_RK2(g,h,r0,v0,t0,tN,N);


%==========
% PLOTTING
%==========
n = linspace(1,N,N);

V1 = zeros(1,length(v1));
% perhaps not the more efficient nor the more elegant way, but we had
% problems using index vectors and 'norm' not returning a vector
for i=1:length(v1)
   V1(i) = norm(v1(:,i)); 
end
I1=2;
while norm(r1(:,I1))>R
    I1 = I1 + 1;
end
V1(I1:length(V1)) = [];
t1 = t;
t1(I1:length(t1)) = [];
plot(t1, V1,'r'), hold on

V2 = zeros(1,length(v2));
for i=1:length(v2)
   V2(i) = norm(v2(:,i)); 
end
I2=2;
while norm(r2(:,I2))>R
    I2 = I2 + 1;
end
V2(I2:length(V2)) = [];
t2 = t;
t2(I2:length(t2)) = [];
plot(t2, V2,'b'), hold on

V3 = zeros(1,length(v3));
for i=1:length(v2)
   V3(i) = norm(v3(:,i)); 
end
I3=2;
while norm(r3(:,I3))>R
    I3 = I3 + 1;
end
V3(I3:length(V3)) = [];
t3 = t;
t3(I3:length(t3)) = [];
plot(t3, V3,'k'), hold on

title('Velocidad en función del tiempo para k1, k2 y k3')
xlabel('Tiempo [s]')
ylabel('Velocidad [m/s]')
legend('k1 = 0.0006','k2 = 0.001','k3 = 0.002')