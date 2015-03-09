% plot Earth
load earth % load image data, X, and colormap, map
sphere; h = findobj('Type','surface');
hemisphere = [ones(257,125),...
              X,...
              ones(257,125)];
set(h,'CData',flipud(hemisphere),'FaceColor','texturemap')
colormap(map)
axis equal
view([90 0])
%%
clear all
clc
close all

R = 6371000;
load earth;
[x,y,z] = sphere(30);

    props.FaceColor= 'texturemap';
    props.EdgeColor = 'none';
    props.FaceLighting = 'phong';
    props.CData = map;

surface(R*x,R*y,R*z,props)
    axis off
view(-180,20)
axis equal