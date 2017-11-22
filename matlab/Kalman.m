function [ ll,store_m ] = Kalman(observations,A,dimension )
 ll = 0;
%  initization
  m = ones(dimension,1);
  P = eye(dimension);
  H = ones(dimension,dimension);
  R = eye(dimension);
  Q = eye(dimension);
  for t = 1 : size(observations,1)
      m_ = A * m;
      P_ = Q+A*P*A';
      mu = H*m_;
      S  = R+H*P_*H'; 
      Vs= chol(S); 
      det_S= prod(diag(Vs))^2; inv_sqrt_S= inv(Vs); iS= inv_sqrt_S*inv_sqrt_S'; 
      K  = P_*H'/S;
      z = observations(t,:)';
      qz_temp =  -0.5*size(z,1)*log(2*pi) - 0.5*log(det_S) - 0.5*dot(z-repmat(mu,[1 size(z,2)]),iS*(z-repmat(mu,[1 size(z,2)])));
      m = repmat(m,[1 size(z,2)]) + K*(z-repmat(mu,[1 size(z,2)]));
      P  = (eye(size(P))-K*H)*P_;
      ll = ll + qz_temp;
      store_m(t,:) = m;
  end
  
end

