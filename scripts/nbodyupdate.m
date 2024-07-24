
function [r_next] = nbodyupdate(m, r_current, r_prev, dt)
% obsolete
% m: Vector of length N containing the particle masses
% r_current: N x 3 array containing the current particle positions
% r_prev: N x 3 array containing the previous particle positions
% dt: number representing time step size

% r_next: N x 3 array containing the computed particle at next time step

% idea: this function outputs an array of the next particle positions,
%       given the current and previous ones. It will use the nbodyaccn 
%       function to determine the accelerations on each particle, then use
%       that to find the future positions, using the formula outlined in
%       the nbody problem explanation.
a = nbodyaccn(m, r_current);
N = size(r_current, 1);
%r_next = zeros(N,3);
%{
for i = 1:N
    % calculate for each dimension
    r_next(i,1) = 2*r_current(i,1) - r_prev(i,1) + a(i,1) * dt^2;
    r_next(i,2) = 2*r_current(i,2) - r_prev(i,2) + a(i,2) * dt^2;
    r_next(i,3) = 2*r_current(i,3) - r_prev(i,3) + a(i,3) * dt^2;
end
%}
r_next = 2.*r_current - r_prev + a.*dt^2;
end