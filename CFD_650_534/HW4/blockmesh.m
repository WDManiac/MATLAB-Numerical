function [dx,dy,mat_ini] = blockmesh(X,Y,gamma,h)
%   Generating rectangular mesh for 2D difference equation solving
%   X, Y ---- total numbers of nodes 
%   gamma ---- ratio of length to height

height = h;
length = gamma*h;

%   dimension of single mesh
dx = length/(X-1);
dy = height/(Y-1);

mat_ini = zeros(Y,X);


