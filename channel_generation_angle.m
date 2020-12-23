function [H,AT,AR,theta_aod,theta_aoa] = channel_generation_angle(N_t,N_r,N_c, N_ray,E_aod, E_aoa, alpha)

%input: the numbers of transmit antennas and receive antennas
%output: the realized channel, codebook for vrf and wrf of OMP method 
%the comments is originally written in Chinese, hopefully I can have time
% to tranfer it to English, but not now.

theta_aoa = repmat(E_aoa,1,N_ray);
aoa = sin(theta_aoa);

theta_aod = repmat(E_aod,1, N_ray);
aod = sin(theta_aod);

signature_t = [0:(N_t-1)]';
signature_t = 1i*pi* signature_t;                          
signature_r = [0:(N_r-1)]';
signature_r = 1i*pi* signature_r;    

H_ray = zeros(N_r, N_t, N_c, N_ray);
H_cl = zeros(N_r, N_t, N_c);


    for i= 1: N_c
        for m = 1: N_ray
            H_ray(:,:,i,m)= alpha(i, m) * exp((aoa(i,m)*signature_r))*exp((aod(i,m)*signature_t))'/sqrt(N_t*N_r); 
        end
    end  
        H_cl = sum(H_ray, 4);    
    
    H(:,:) = sqrt(N_t*N_r/N_c/N_ray)*sum(H_cl,3);   
    
    aod = aod(:).';
    aoa = aoa(:).';
    A = kron(aod,signature_t);
    AT = 1/sqrt(N_t)*exp(A);
    A = kron(aoa,signature_r);
    AR = 1/sqrt(N_r)*exp(A);
