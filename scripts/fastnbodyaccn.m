
function [result] = fastnbodyaccn(m, stars, cores)
% m: Vector of length ncore containing the core masses
% stars: N x 3 array containing the particle positions
% cores: ncores x 3 array containing the core positions
% a: (N+ncore) x 3 array containing the computed particle and core 
%     accelerations

% idea: this function outputs an array of the corresponding acceleration
%       components for each dimension, for each particle. It is more
%       efficient because other stars do not effect the acceleration
%       calculations. 

% number of particles
N = size(stars, 1);
ncore = size(cores,1);

% acceleration vector for stars and cores
[a] = zeros(N,3);
core_acc = zeros(ncore,3);

% loop over all cores
for i = 1:ncore
    % calculate the separation
    corestack = repmat(cores(i, :), N, 1);
    deltas = corestack - stars;
    r_ij_3 = (deltas(:, 1).^2 + deltas(:, 2).^2 + deltas(:, 3).^2).^(3/2);

    % make an array the same sie as the accel vector, but with just the
    % position of a core
    masses = repmat(m(i), N, 1);
    % update acceleration vector through array operations
    a = a + masses .* deltas ./ r_ij_3;

    % loop over all other cores to calculate the accel. due to other cores
    for j = 1 : ncore
        % excluding itself
        if i == j
            % don't put break here, as it will stop the whole for loop
            % just want to skip one instance. 
            %disp('do nothing')
        else 
            % calculate the separation
            delta_x = cores(j,1) - cores(i,1);
            delta_y = cores(j,2) - cores(i,2);
            delta_z = cores(j,3) - cores(i,3);

            r_ij_3 = ( (delta_x)^2 + (delta_y)^2 + (delta_z)^2 )^(3/2);
            % update acceleration vector through the sum
            core_acc(i,1) = core_acc(i,1) + m(j) * delta_x / r_ij_3;
            core_acc(i,2) = core_acc(i,2) + m(j) * delta_y / r_ij_3;
            core_acc(i,3) = core_acc(i,3) + m(j) * delta_z / r_ij_3;
        end
    end
end
result = cat(1, core_acc, a);
end