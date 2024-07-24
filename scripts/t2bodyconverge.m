
% script to test convergence of our FDA and equations of motion
% uses twobodysim.m to run the simulations

% Mass of particles
m = [1 0.5];

% particle separation
separation = 16;

% max time
tmax = 60;

% number of particles
N = 2;

% run for levels of interest
[r6, x6, t6] = twobodysim(m, tmax, separation, N, 6);
[r7, x7, t7] = twobodysim(m, tmax, separation, N, 7);
[r8, x8, t8] = twobodysim(m, tmax, separation, N, 8);
[r9, x9, t9] = twobodysim(m, tmax, separation, N, 9);

% prepare subplots
fig = tiledlayout(3,1);
fontsize = 12;

% plot the x values
nexttile
hold on; 
plot(t6, x6, 'r-.o');
plot(t7, x7, 'g-.+'); 
plot(t8, x8, 'b-.*');
plot(t9, x9, '--');
xlabel('Time step','FontSize',fontsize)
ylabel('x position','FontSize',fontsize)
title('x-values of 1 Particle (Levels 6, 7, 8, 9)','FontSize',fontsize)

% downsample
x7 = x7(1:2:end);
x8 = x8(1:4:end);
x9 = x9(1:8:end);

% find difference
dx67 = x6 - x7;
dx78 = x7 - x8;
dx89 = x8 - x9;

% plot the differences
nexttile
hold on; 
plot(t6, dx67, 'r-.o'); 
plot(t6, dx78, 'g-.+');
plot(t6, dx89, '--');
xlabel('Time step','FontSize',fontsize)
ylabel('x position','FontSize',fontsize)
title('Differences of Downsampled x-values','FontSize',fontsize)

% scale them and see if they are coincident
dx78 = 4 * dx78;
dx89 = 16 * dx89;

% plot again
nexttile
hold on; 
plot(t6, dx67, 'r-.o'); 
plot(t6, dx78, 'g-.+');
plot(t6, dx89, '--');
xlabel('Time step','FontSize',fontsize)
ylabel('x position','FontSize',fontsize)
title('Scaled Differences of Downsampled x-values','FontSize',fontsize);
%savefig(fig, 'convergence_testing.png');