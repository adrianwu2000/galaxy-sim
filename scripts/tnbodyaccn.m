
% test case for nbodyaccn in the scenario where 2 bodies orbit about their
% center of mass
m = [1 0.5];
% distance between the two masses
separation = 16;
% Calculating time steps
tmax = 350;
level = 8;
dt = tmax* 2^(-level);
nt = 2^level + 1;
t = linspace(0.0, tmax, nt);
% number of particles
N = 2;
% array to store all the data
r = zeros(N, 3, nt);
% initial conditions for particle 1 at first time step
r1 = m(2) * separation / (m(1) + m(2));
r2 = m(1) * separation / (m(1) + m(2));
r(1, 1, 1) = r1;
r(1, 2, 1) = 0;
r(1, 3, 1) = 0;
% initial conditions for particle 2 at first time step
r(2, 1, 1) = - r2;
r(2, 2, 1) = 0;
r(2, 3, 1) = 0;
% initial velocities
v1 = sqrt(m(2)*r1) / separation;
v2 = sqrt(m(1)*r2) / separation;
v = [0 v1 0; 0 -v2 0];
% acceleration based on initial values
acc = nbodyaccn(m, r(:, :, 1));
% use taylor expansion to get the values at the second time step
r(1, 1, 2) = r(1, 1, 1) + dt * v(1, 1) + 0.5 * dt^2 * acc(1, 1);
r(1, 2, 2) = r(1, 2, 1) + dt * v(1, 2) + 0.5 * dt^2 * acc(1, 2);
r(1, 3, 2) = r(1, 3, 1) + dt * v(1, 3) + 0.5 * dt^2 * acc(1, 3);

r(2, 1, 2) = r(2, 1, 1) + dt * v(2, 1) + 0.5 * dt^2 * acc(2, 1);
r(2, 2, 2) = r(2, 2, 1) + dt * v(2, 2) + 0.5 * dt^2 * acc(2, 2);
r(2, 3, 2) = r(2, 3, 1) + dt * v(2, 3) + 0.5 * dt^2 * acc(2, 3);

% prepare animation
% enable/disable plotting
plotenable = 1;

% speed of animation
pausesecs = 0.0;

% Plot attributes defining the appearance  of the planet
% Particle has a (marker) size of 15 ...

particlesize = 10;
particlecolor = 'r';
particlemarker = '.';

% Size of window for animation
mlim = 15;
% border width
dlim = 0.2*mlim;

% Set avienable to a non-zero value to make an AVI movie.
avienable = 1;

% If plotting is disabled, ensure that AVI generation
% is as well
if ~plotenable
   avienable = 0;
end

% Presumed AVI playback rate in frames per second.
aviframerate = 25;
avifilename = 'two_body_orbit.avi';

% Initialize an avi object.
if avienable
   aviobj = VideoWriter(avifilename);
   open(aviobj);
end

% iterate over all the time steps and update the solution vector r
% we start at n+1 = 3, until the max number of steps, nt, is reached
clf;
for step = 2:nt-1
    if plotenable
      hold on;

      % Define plotting area
      axis square;
      box on;
      xlim([-mlim-dlim, mlim + dlim]);
      ylim([-mlim-dlim, mlim + dlim]);

      % title. 
      titlestr = sprintf('Step: %d', step);
      title(titlestr, 'FontSize', 16, 'FontWeight', 'bold', ...
         'Color', [0.25, 0.42, 0.31]);

      % Draw the particles
      for k = 1:N
          plot(r(k,1,step), r(k, 2, step), 'Marker', particlemarker, 'MarkerSize', particlesize, ...
         'MarkerEdgeColor', particlecolor, 'MarkerFaceColor', particlecolor)
      end

      drawnow;
    
      % Record video frame if AVI recording is enabled and record 
      % multiple copies of the first frame
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
    % find current and previous slices
    r_current = r(:, :, step);
    r_prev = r(:, :, step-1);
    % Get the acceleration 
    a = nbodyaccn(m, r_current);
    % Get next positions
    r_next = 2.*r_current - r_prev + a.*dt^2;
    r(:, :, step+1) = r_next;
end

% write a diagnostic message that a movie file was created.
if avienable
   close(aviobj);
   fprintf('Created video file: %s\n', avifilename);
end