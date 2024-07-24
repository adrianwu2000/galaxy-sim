
% test case for altered version of nbodyaccn.m with 2 galaxiesnbodyaccn
% make results repeatable
rng("default");

% two center of masses and many stars around it
ncore = 2;
mc1 = 8; % mass of core
mc2 = 8;

% time step stuff
tmax = 400;
level = 8;
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
core1 = [-40, 30, 0];
core2 = [0, -35, 0];
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
kickx = 0.5;
kicky = 0.2;
% core1
v(ncore + 1:N + ncore, 1) = v1 .* cos(theta1 + pi/2) + kickx;
v(ncore + 1:N + ncore, 2) = v1 .* sin(theta1 + pi/2) + kicky;
% core2
v(1 + N + ncore:end, 1) = v2 .* cos(theta2 - pi/2) -0 ;
v(1 + N + ncore:end, 2) = v2 .* sin(theta2 - pi/2) -0 ;

% core velocities
v(1, 1) = kickx;
v(1, 2) = kicky;
v(2, 1) = -0;
v(2, 2) = -0;

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

%-----------------------------------------------------------
% Set plotenable to non-zero/zero to enable/disable plotting.
%-----------------------------------------------------------
plotenable = 1;
%-----------------------------------------------------------
% Parameter to control speed of animation.  Script execution
% will pause for pausesecs each time a new frame is drawn.
% 
% Setting this parameter to a "largish" value, say 0.1
% (seconds), will produce a slow-motion effect.
% 
% Set it to 0 for maximum animation speed.
%-----------------------------------------------------------
pausesecs = 0.0;
%-----------------------------------------------------------
% Plot attributes defining the appearance  of the planet.
%-----------------------------------------------------------

% Ball has a (marker) size of
ballsize = 8;
% ... it's red ...
ballcolor = 'r';
% ... and it's plotted as a circle.
ballmarker = '.';
%-----------------------------------------------------------
mlim = 100;
dlim = 0.2*mlim;
%-----------------------------------------------------------
% Set avienable to a non-zero value to make an AVI movie.
%-----------------------------------------------------------
avienable = 1;

% If plotting is disabled, ensure that AVI generation
% is as well
if ~plotenable
   avienable = 0;
end

% Name of avi file.
%avifilename = 'bounce.avi';

% Presumed AVI playback rate in frames per second.
aviframerate = 25;
avifilename = 'two_galaxy_sim_tester.avi';

%-----------------------------------------------------------
% If AVI creation is enabled, then initialize an avi object.
%-----------------------------------------------------------
if avienable
   aviobj = VideoWriter(avifilename);
   open(aviobj);
end

%-----------------------------------------------------------
% SIMULATE
%-----------------------------------------------------------

% iterate over all the time steps and update the solution vector r
% note we start at n+1 = 3, until the max number of steps, nt, is reached
for step = 2:nt-1
    if plotenable
      % Clear figure
      clf;

      % Don't erase figure after each plot command.
      hold on;

      % Define plotting area, using a 1:1 aspect ratio for the 
      % plotted region, boxed axes and a 15%-width "border" around 
      % the unit square.
      axis square;
      box on;
      xlim([-mlim-dlim, mlim + dlim]);
      ylim([-mlim-dlim, mlim + dlim]);

      % Make and display title. 
      titlestr = sprintf('Step: %d | Stars per Core: %d', step, N);
      %titlestr = sprintf('Step: %d X: %.1f Y: %.1f', step, ...
      %    r(1,1,step), r(1,2,step));
      title(titlestr, 'FontSize', 16, 'FontWeight', 'bold', ...
         'Color', [0.25, 0.42, 0.31]);

      % Draw the ball.
      %disp('plot begin')
      plot(r(1,1,step), r(1, 2, step), 'Marker', 'o', 'MarkerSize', 15, ...
         'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');
      plot(r(2,1,step), r(2, 2, step), 'Marker', 'o', 'MarkerSize', 15, ...
         'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');
      plot(r(ncore+1:ncore+N, 1, step), r(ncore+1:ncore+N, 2, step), ballmarker, 'MarkerSize', ballsize, ...
         'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');
      plot(r(ncore+N+1:end, 1, step), r(ncore+1+N:end, 2, step), ballmarker, 'MarkerSize', ballsize, ...
         'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');      
        %{
      for k = ncore+1:size(r,1)
          ballcolor = 'r';
          if k > N + ncore
              ballcolor = 'b';
          end
          plot(r(k,1,step), r(k, 2, step), 'Marker', ballmarker, 'MarkerSize', ballsize, ...
         'MarkerEdgeColor', ballcolor, 'MarkerFaceColor', ballcolor)
      end
        %}
      %disp('plot end')
      % For
      % Force update of figure window.
      drawnow;
    
      % Record video frame if AVI recording is enabled and record 
      % multiple copies of the first frame.  Here we record 5 seconds
      % worth which will allow the viewer a bit of time to process 
      % the initial scene before the animation starts.
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
    %disp('dynamics start')
    % dynamics
    r_current = r(:, :, step);
    r_prev = r(:, :, step-1);
    r_next = fastnbodyupdate(N, ncore, m, r_current, r_prev, dt);
    r(:, :, step+1) = r_next;
    %disp('dynamics end')
end
% If we're making a video, finalize the recording and 
% write a diagnostic message that a movie file was created.

if avienable
   close(aviobj);
   fprintf('Created video file: %s\n', avifilename);
end