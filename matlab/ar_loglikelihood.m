function [ ll ] = ar_loglikelihood( xparticles,z,R )
M = size(xparticles,2);
% e_sq= sum( (diag(1./diag(R))*(repmat(z',[1 M])- xparticles)).^2 );
% ll =  -e_sq/2 -  log(prod(diag(R))) - size(xparticles,1)* 0.9189385  ;
for i = 1 : M
    e_sq =  (z - xparticles(:,i))*inv(R)* (z - xparticles(:,i));
    ll(i,1) = -e_sq/2 -  log(prod(diag(R)))- size(xparticles,1)* 0.9189385 ;
end
end

