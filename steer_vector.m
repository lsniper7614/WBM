function [steering_vector] = steer_vector(N, theta)
%steer_vector : Generate the antenna array vector 
%   N   £º The number of antenna in ULA
% theta £º  The AoA or AoD
signature = 1j * pi * (0:N-1)';
steering_vector = exp(signature .* sin(theta.'));
end

