
function [r_next] = fastnbodyupdate(N, ncore, m, r_current, r_prev, dt)
% obsolete

% N: number of particles (stars) per core
% ncore: number of cores
% m: Vector of length ncore containing the core masses
% r_current: (N+ncore) x 3 array containing the current particle positions
% r_prev: (N+ncore) x 3 array containing the previous particle positions
% dt: number representing time step size

% r_next: (N+ncore) x 3 array containing the computed particle at next time step

% idea: this function outputs an array of the next particle positions,
%       given the current and previous ones. It will use the fastnbodyaccn 
%       function to determine the accelerations on each particle, then use
%       that to find the future positions, using the formula outlined in
%       the nbody problem explanation. The stars won't affect the
%       acceleration and position changes. 
%
%       it also assumes that the first ncore rows of r_current are the
%       positions of the cores, and that the rest are the position of the
%       stars.

% find star and core positions from input arrays
star_current = r_current(ncore+1:end, :);
%disp("size of star current is")
%disp(size(star_current))

core_current = r_current(1:ncore, :);
%disp("size of core current is")
%disp(size(core_current))
%star_prev = r_prev(ncore+1:end);
%core_prev = r_prev(1:ncore);

% find current acceleration values
a = fastnbodyaccn(m, star_current, core_current);

% calculate next set of positions based on FDA
%r_next = zeros(ncore*N + ncore,3);
%for i = 1:ncore*N+ncore
%    % calculate for each dimension
%    r_next(i,1) = 2*r_current(i,1) - r_prev(i,1) + a(i,1) * dt^2;
%    r_next(i,2) = 2*r_current(i,2) - r_prev(i,2) + a(i,2) * dt^2;
%    r_next(i,3) = 2*r_current(i,3) - r_prev(i,3) + a(i,3) * dt^2;
%end
r_next = 2.*r_current - r_prev + a.*dt^2;