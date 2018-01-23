function [x, P, VV] = kalman_f(y, A, C, Q, R, init_x, init_P)
% Kalman filter 
%2018.1
% [x, P, VV] = kalman_f(y, A, C, Q, R, init_x, init_P)
%
% INPUTS:
% x(t)=A*x(t-1)+Q;
% y(t)=C*x(t)+R
% init_x - the initial state (column) vector 
% init_P- the initial state covariance 

% OUTPUTS 
% x(:,t) = E[X(:,t) | y(:,1:t)]
% P(:,:,t) = Cov[X(:,t) | y(:,1:t)]
% VV(:,:,t) = Cov[X(:,t), X(:,t-1) | y(:,1:t)] t >= 2


[os T] = size(y);
ss = size(A,1); % size of state space

model = ones(1,T);

x = zeros(ss, T);
P = zeros(ss, ss, T);
VV = zeros(ss, ss, T);


for t=1:T
  m = model(t);
  if t==1
    %prevx = init_x(:,m);
    %prevV = init_V(:,:,m);
    prevx = init_x;
    prevP = init_P;
    initial = 1;
  else
    prevx = x(:,t-1);
    prevP = P(:,:,t-1);
    initial = 0;
  end


    xpred = prevx;
  
if initial
  Vpred = prevP;
 else
  Vpred = A(:,:,m)*prevP*A(:,:,m)' + Q(:,:,m);
end

e = y(:,t) - C(:,:,m)*xpred; % error (innovation)
n = length(e);
ss = length(A(:,:,m));
S = C(:,:,m)*Vpred*C(:,:,m)' + R(:,:,m);
Sinv = inv(S);
ss = length(prevP);

K = Vpred*C(:,:,m)'*Sinv; % Kalman gain matrix
% If there is no observation vector, set K = zeros(ss).
x(:,t) = xpred + K*e;
P(:,:,t)= (eye(ss) - K*C(:,:,m))*Vpred;
VV(:,:,t) = (eye(ss) - K*C(:,:,m))*A(:,:,m)*prevP;

  
  end

end





