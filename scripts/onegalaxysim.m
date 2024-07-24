
% test case for fastnbodyaccn.m

% make results repeatable
rng("default");

% two center of masses and many stars around it
ncore = 1;
mc1 = 9; % mass of core

% time step stuff
tmax = 100;
level = 7;
dt = tmax* 2^(-level);
nt = 2^level + 1;
t = linspace(0.0, tmax, nt);

% number of particles per galaxy
N = 5000;

% mass of all stars and core, mass of core is always first, mass of stars
% is zero
m = mc1;

% array to store all the data
r = zeros(N + ncore, 3, nt);

% min and max radii about the core
rmin = 2;
rmax = 5;
% core initial positions
core1 = [10, 10, 0];
% intial positions (assumes z=0)
[init1, theta1] = randcirclepts(rmin, rmax, N, core1(1), core1(2));
r(:, :, 1) = cat(1, core1, init1);

%%%%%%%%%%%%%%%%%%%%%%%

% initial velocities
v = zeros(N + ncore, 3);

% find intial velocity based on separation from core
% skip over core, because it is stationary

% assume 2D
separation1 = sqrt((r(:, 1) - core1(1)).^2 + ...
    (r(:, 2) - core1(2)).^2);

% inital velocity for circular orbit
v0 = sqrt(mc1 ./ separation1);

% find angle of each star
theta1 = cat(1, 0, theta1);

% initial kick for translation
kickx = 0.2;
kicky = 0.2;

% core1 inital velocities for stars
% change sign of pi/2 to change rotation direction
v(:, 1) = v0 .* cos(theta1 - pi/2) - kickx;
v(:, 2) = v0 .* sin(theta1 - pi/2) - kicky;

% core velocities
v(1, 1) = -kickx;
v(1, 2) = -kicky;

%%%%%%%%%%%%%%%%%%%%%%%%
% acceleration based on initial values
acc = fastnbodyaccn(m, r(2:end, :, 1), r(1, :, 1));

% use taylor expansion to get the values at the second time step
r(:, :, 2) = r(:, :, 1) + dt .* v(:, :) + 0.5 * dt^2 .* acc(:, :);

% prepare animation

% enable plotting 
plotenable = 1;
% speed of plot
pausesecs = 0.0;

% Plot variables
particlesize = 8;
particlecolor = 'r';
particlemarker = '.';

% size of animation window and border
mlim = 20;
dlim = 0.2*mlim;

% enable animation
avienable = 1;

% If plotting is disabled, ensure that AVI generation
% is as well
if ~plotenable
   avienable = 0;
end

% Presumed AVI playback rate in frames per second.
aviframerate = 25;
avifilename = 'one_galaxy_sim.avi';

%initialize an avi object.
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
      hold on;

      % Define plotting area
      axis square;
      box on;
      xlim([-mlim-dlim, mlim + dlim]);
      ylim([-mlim-dlim, mlim + dlim]);

      % Make and display title. 
      titlestr = sprintf('Step: %d', step);
      title(titlestr, 'FontSize', 16, 'FontWeight', 'bold', ...
         'Color', [0.25, 0.42, 0.31]);

      % Plot the particles
      plot(r(1,1,step), r(1, 2, step), 'Marker', particlemarker, 'MarkerSize', 10, ...
         'MarkerEdgeColor', 'c', 'MarkerFaceColor', 'c')
      plot(r(2:end, 1, step), r(2:end, 2, step), particlemarker, 'MarkerSize', particlesize, ...
         'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');

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
    % get current and previous positions
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