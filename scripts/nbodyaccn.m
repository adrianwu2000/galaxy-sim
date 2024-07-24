
function [a] = nbodyaccn(m, r)
% m: Vector of length N containing the particle masses
% r: N x 3 array containing the particle positions
% a: N x 3 array containing the computed particle accelerations

% idea: this function outputs an array of the corresponding acceleration
%       components for each dimension, for each particle.

N = size(r, 1);
[a] = zeros(N,3);
% for each particle
for i = 1: N
    % loop over all other particles
    for j = 1:N
        % excluding itself
        if i == j
            % don't put break here, as it will stop the whole for loop
            % just want to skip one instance. 
            %disp('do nothing')
        else 
            % calculate the separation
            delta_x = r(j,1) - r(i,1);
            delta_y = r(j,2) - r(i,2);
            delta_z = r(j,3) - r(i,3);

            r_ij_3 = ( (delta_x)^2 + (delta_y)^2 + (delta_z)^2 )^(3/2);
            % update acceleration vector through the sum
            a(i,1) = a(i,1) + m(j) * delta_x / r_ij_3;
            a(i,2) = a(i,2) + m(j) * delta_y / r_ij_3;
            a(i,3) = a(i,3) + m(j) * delta_z / r_ij_3;
        end
    end
end
end