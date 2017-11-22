function [index] = stratified_resampling(wn,N2,st)

% wn: weight vector
% N2: number of particles to sample

if (nargin==3)
    rand('state',st);
end;

wn = wn(:);
N = size(wn,1);

wn = wn/sum(wn);


g = rand/N2;

r = [g:1/N2:1+g]';

jj = 1;

x2 = cumsum(wn);

for ii = 1:N
    while x2(ii)>r(jj)
        index(jj) = ii;
        jj = jj+1;
    end;
end;

index = index(:);