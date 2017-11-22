function [ ll ] = ar_loglikelihood( xparticles,z,R )
M = size(xparticles,2);
e_sq= sum( (diag(1./diag(R))*(repmat(z',[1 M])- xparticles)).^2 );
ll =  -e_sq/2 - log(2*pi*prod(diag(R))) ;
end

