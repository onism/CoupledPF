function [ particles ] = ar_transition( particles,A,randomness)
%  [no,dim] = size(particles);
%  particles = A * particles + gen_gms(1,zeros(0,dim),eye(dim),no);
 particles = A * particles + randomness;
end

