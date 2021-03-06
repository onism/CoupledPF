clear  
close all
clc
dbstop if error
rand('seed',777)
dimension = 1;
alpha_star = 0.95;
A_star = createA(alpha_star,dimension);
x_t = gen_gms(1,zeros(dimension,1),eye(dimension),1);
datalength = 100;
for t = 1 : datalength
    x_t = A_star * x_t + gen_gms(1,zeros(dimension,1),eye(dimension),1);
    observations(t,:) = x_t + gen_gms(1,zeros(dimension,1),eye(dimension),1);
    store_x(t,:) = x_t;
end


% now we estimat the loglihood with various methods
nparticles = 2^7;
nrep = 5;
nrhos = 20;
rhos = linspace(0.3,0.4,20);
%  estimate the log-likelihood using independent bootstrap particle filters
% parfor i = 1 : nrep
for i = 1 : nrep
    count = 1;
    for irho = rhos
        theta = irho;
        for t = 1 : datalength+1
            randomness(t,:) =   gen_gms(1,zeros(dimension,1),eye(dimension),nparticles);
        end
        
        ll = particle_filter_storeall(nparticles,theta,observations,randomness,dimension);
        pfll(i,count) = ll;
        count = count + 1;
    end
end


%   Common random numbers

for i = 1 : nrep
    count = 1;
    for t = 1 : datalength+1
            randomness(t,:) =   gen_gms(1,zeros(dimension,1),eye(dimension),nparticles);
    end
%     randomness =   gen_gms(1,zeros(dimension,1),eye(dimension),nparticles); 
    for irho = rhos
        theta = irho;
        ll = particle_filter_storeall(nparticles,theta,observations,randomness,dimension);
        common_pfll(i,count) = ll;
        count = count + 1;
    end
end
 

% evaluate the log-likelihood exactly using the Kalman filter
 
count = 1;
for irho = rhos
   theta = irho;
   A = createA(irho,dimension);
   ll = Kalman(observations,A,dimension );
   kfll(count,1) = ll;
   count = count + 1;
end
 
 
[ll,store_m] = Kalman(observations,A,dimension );
 

figure
hold on
for i = 1 : nrep
    plot(pfll(i,:));
end
% plot(kfll)


figure
hold on
for i = 1 : nrep
    plot(common_pfll(i,:));
end
% plot(kfll)


% figure
% hold on
% plot(store_m(1:10,2))
% plot(store_x(1:10,2))
