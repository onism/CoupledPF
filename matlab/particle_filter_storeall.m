function [ ll ] = particle_filter_storeall( nparticles,theta,observations, randomness,dimensional)
 datalength = size(observations,1);
 A = createA(theta,dimensional);
%  init particles
 xparticles = randomness;
 normweights = 1/ nparticles * ones(1,nparticles);
 ll = 0;
 R = eye(dimensional);
 
%  begin filter
    for t =1 : datalength
%         resampling
        ancestors = stratified_resampling(normweights,nparticles);
        xparticles = xparticles(:,ancestors);
%         prediction
        xparticles = ar_transition(xparticles,A,randomness);
%         compute log likelihood
        logw = ar_loglikelihood( xparticles,observations(t,:),R );
        maxlw = max(logw);
        w = exp(logw - maxlw);
        ll = ll + maxlw +  log(mean(w));
        normweights = w / sum(w);
    end
end

