function [Wrf] = beam_alignment(recv_sig, recv_antenna)
% beam_alignment 

scan_theta = -90: 1 :90;
scan_vec = steer_vector(recv_antenna, scan_theta.');
scan_y = scan_vec' * recv_sig;
scan_y_energy = mean(abs(scan_y), 2);
[~, Wrf_index] = max(scan_y_energy);
Wrf = scan_vec(:, Wrf_index);
end

