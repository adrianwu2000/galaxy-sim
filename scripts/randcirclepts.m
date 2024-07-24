
function [pos, theta] = randcirclepts(rmin, rmax, n, x0, y0)
% rmin: minimum radius for this circle (must be greater than 0)
% rmax: maximum radius for this circle
%$ n: how many random points to create
% x0: position of center for this circle
% y0: position of center for this circle
% r: an array of x,y,z coordinates of the positions of the random points

% this function generates the array of random star positions for the galaxy
% generation over a disk or donut shape. It only distributes over a 2D
% plane (x,y), but will output a third dimension for the sake of the
% project (z).

% random angles 
theta = 2*pi*rand(n,1);

% random radii
r = (rmax-rmin)*rand(n,1) + rmin;

% combine them to get the random positions of stars
x = x0 + r.*cos(theta);
y = y0 + r.*sin(theta);

% optional plot to see how it looks
%plot(x,y, '.');

% combine the x, y positions
z = zeros(n,1);
pos = cat(2, x, y, z);

end