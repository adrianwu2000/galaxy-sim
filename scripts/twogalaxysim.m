
% 2 galaxies collide

% make results repeatable
rng("default");

% two center of masses and many stars around it
ncore = 2;
mc1 = 8; % mass of core
mc2 = 8;

% time step stuff
tmax = 400;
level = 11;
dt = tmax* 2^(-level);
nt = 2^level + 1;
t = linspace(0.0, tmax, nt);

% number of particles per galaxy
N = 10000;

% mass of all stars and core, mass of core is always first, mass of stars
% is zero
m = cat(2, mc1, mc2, zeros(1,2*N));

% array to store all the data
r = zeros(2*N + ncore, 3, nt);

% min and max radii about the core
rmin = 4;
rmax = 30;

% core initial positions
core1 = [-40, 10, 0];
core2 = [40, -10, 0];

% intial positions (assumes z=0)
[init1, theta1] = randcirclepts(rmin, rmax, N, core1(1), core1(2));
[init2, theta2] = randcirclepts(rmin, rmax, N, core2(1), core2(2));
r(:, :, 1) = cat(1, core1, core2, init1, init2);

%%%%%%%%%%%%%%%%%%%%%%%

% initial velocities
v = zeros(2*N + ncore, 3);

% find intial velocity based on separation from core
% skip over core, because it is stationary

% assume 2D
separation1 = sqrt((r(ncore + 1:N + ncore, 1) - core1(1)).^2 + ...
    (r(ncore + 1:N + ncore, 2) - core1(2)).^2);
separation2 = sqrt((r(1 + N + ncore:end, 1) - core2(1)).^2 + ...
    (r(1 + N + ncore:end, 2) - core2(2)).^2);

% inital velocity for circular orbit
v1 = sqrt(mc1 ./ separation1);
v2 = sqrt(mc2 ./ separation2);
v0 = cat(1, v1, v2);

% initial kick for translation
kickx = 0.2;
kicky = 0.3;

% core1 stars initial velocities
% (change sign of pi/2 to change orbit direction)
v(ncore + 1:N + ncore, 1) = v1 .* cos(theta1 + pi/2) + kickx;
v(ncore + 1:N + ncore, 2) = v1 .* sin(theta1 + pi/2) + kicky;

% core2 stars initial velocities
% (change sign of pi/2 to change orbit direction)
v(1 + N + ncore:end, 1) = v2 .* cos(theta2 - pi/2) - kickx;
v(1 + N + ncore:end, 2) = v2 .* sin(theta2 - pi/2) - kicky;

% core velocities
v(1, 1) = kickx;
v(1, 2) = kicky;
v(2, 1) = -kickx;
v(2, 2) = -kicky;

%%%%%%%%%%%%%%%%%%%%%%%%
% acceleration based on initial values
acc = fastnbodyaccn(m, r(ncore+1:end, :, 1), r(1:ncore, :, 1));

% use taylor expansion to get the values at the second time step
for i = 1:size(r,1)
    r(i, 1, 2) = r(i, 1, 1) + dt * v(i, 1) + 0.5 * dt^2 * acc(i, 1);
    r(i, 2, 2) = r(i, 2, 1) + dt * v(i, 2) + 0.5 * dt^2 * acc(i, 2);
    r(i, 3, 2) = r(i, 3, 1) + dt * v(i, 3) + 0.5 * dt^2 * acc(i, 3);
end

% prepare animation

% enable plotting
plotenable = 1;
% pause time
pausesecs = 0.0;

particlesize = 8;
particlecolor = 'r';
ballmarker = '.';

% window for animation and border
mlim = 100;
dlim = 0.2*mlim;

% make a movie
avienable = 1;

% If plotting is disabled, ensure that AVI generation
% is as well
if ~plotenable
   avienable = 0;
end

% Presumed AVI playback rate in frames per second.
aviframerate = 25;
avifilename = 'two_galaxy_sim.avi';

if avienable
   aviobj = VideoWriter(avifilename);
   open(aviobj);
end

% iterate over all the time steps and update the solution vector r
% note we start at n+1 = 3, until the max number of steps, nt, is reached
for step = 2:nt-1
    if plotenable
      % Clear figure
      clf;

      % Don't erase figure after each plot command.
      hold on;

      % Define plotting area
      axis square;
      box on;
      xlim([-mlim-dlim, mlim + dlim]);
      ylim([-mlim-dlim, mlim + dlim]);

      % Make and display title. 
      titlestr = sprintf('Step: %d | Stars per Core: %d', step, N);
      title(titlestr, 'FontSize', 16, 'FontWeight', 'bold', ...
         'Color', [0.25, 0.42, 0.31]);

      % plot stars and cores
      plot(r(1,1,step), r(1, 2, step), 'Marker', 'o', 'MarkerSize', 15, ...
         'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');
      plot(r(2,1,step), r(2, 2, step), 'Marker', 'o', 'MarkerSize', 15, ...
         'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');
      plot(r(ncore+1:ncore+N, 1, step), r(ncore+1:ncore+N, 2, step), ballmarker, 'MarkerSize', particlesize, ...
         'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');
      plot(r(ncore+N+1:end, 1, step), r(ncore+1+N:end, 2, step), ballmarker, 'MarkerSize', particlesize, ...
         'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');      
        
      % Force update of figure window.
      drawnow;
    
      if avienable
         if t == 0
            framecount = 5 * aviframerate ;
         else
            framecount = 1;
         end
         for iframe = 1 : framecount
            writeVideo(aviobj, getframe(gcf));
         end
      end
      % Pause execution to control interactive visualization speed.
      pause(pausesecs);
    end

    % dynamics

    % find current and previous positions
    r_current = r(:, :, step);
    r_prev = r(:, :, step-1);

    % split into stars vs. cores
    star_current = r_current(ncore+1:end, :);
    core_current = r_current(1:ncore, :);
    
    % find current acceleration values
    a = fastnbodyaccn(m, star_current, core_current);
    
    % calculate next set of positions based on FDA
    r_next = 2.*r_current - r_prev + a.*dt^2;
    r(:, :, step+1) = r_next;

end

if avienable
   close(aviobj);
   fprintf('Created video file: %s\n', avifilename);
end