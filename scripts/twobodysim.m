
function [r, x, t] = twobodysim(m, tmax, separation, N, level)
% N: number of particles (probably 2)
% m: length N array that has the masses of each particle
% tmax: maximum time
% separation: distance between the particles, assumes they are on the xaxis
% l: the level which defines how many time steps to use

% r: a Nx3x "number of time steps" array that contains the positions of all
%   the particles at all time steps
% x: a size 1x "number of time steps" array that has all the x positions of
%   one particle of choice (arbitrary)

% time steps
dt = tmax* 2^(-level);
nt = 2^level + 1;

% array to store all the positions
r = zeros(N, 3, nt);
x = zeros(1, nt);
t = linspace(0.0, tmax, nt);

% initial conditions for particle 1 at first time step
r1 = m(2) * separation / (m(1) + m(2));
r2 = m(1) * separation / (m(1) + m(2));
r(1, 1, 1) = r1;
r(1, 2, 1) = 0;
r(1, 3, 1) = 0;
x(1,1) = r1;

% initial conditions for particle 2 at first time step
r(2, 1, 1) = - r2;
r(2, 2, 1) = 0;
r(2, 3, 1) = 0;
%x(2, 1) = -r2;

% initial velocities for orbital motion
v1 = sqrt(m(2)*r1) / separation;
v2 = sqrt(m(1)*r2) / separation;
v = [0 v1 0; 0 -v2 0];

% initial acceleration based on initial condition
acc = nbodyaccn(m, r(:, :, 1));

% use taylor expansion to get the values at the second time step
r(1, 1, 2) = r(1, 1, 1) + dt * v(1, 1) + 0.5 * dt^2 * acc(1, 1);
r(1, 2, 2) = r(1, 2, 1) + dt * v(1, 2) + 0.5 * dt^2 * acc(1, 2);
r(1, 3, 2) = r(1, 3, 1) + dt * v(1, 3) + 0.5 * dt^2 * acc(1, 3);
x(1,2) = r(1, 1, 2);

r(2, 1, 2) = r(2, 1, 1) + dt * v(2, 1) + 0.5 * dt^2 * acc(2, 1);
r(2, 2, 2) = r(2, 2, 1) + dt * v(2, 2) + 0.5 * dt^2 * acc(2, 2);
r(2, 3, 2) = r(2, 3, 1) + dt * v(2, 3) + 0.5 * dt^2 * acc(2, 3);
%x(2,2) = r(2, 1, 2);

% iterate over all the time steps and update the solution vector r
% we start at n+1 = 3, until the max number of steps is reached
for step = 2:nt-1
    % dynamics
    % find current and previous position slices
    r_current = r(:, :, step);
    r_prev = r(:, :, step-1);
    % get acceleration
    a = nbodyaccn(m, r_current);
    % get the next positions
    r_next = 2.*r_current - r_prev + a.*dt^2;
    % add to our solution vectors
    r(:, :, step+1) = r_next;
    x(1, step+1) = r_next(1,1);
end
end