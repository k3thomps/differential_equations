function [cost,grad] = costFunction(in,L,n,lambda_f,lambda_s,beta_,phase,R)

% initialization

dx = (2*L)/n;

f = in(1:n+1);
s = in(n+2:end);

% cost function

cost = 0;
for i = 2:n,
 cost = cost + (1/(2*dx))*(f(i) - f(i-1))^2 + (1/(2*dx))*(s(i) - s(i-1))^2 + lambda_f*(f(i)^2 - 1)^2*dx + lambda_s*(s(i)^2 - 2)*s(i)^2*dx + beta_*f(i)^2*s(i)^2*dx + (dx/2)*(phase/R)^2*s(i)^2;
endfor;

% gradient

delta = zeros(n+1,2);

for i=2:n,
 delta(i,1) = (-1/(dx^2))*(f(i+1) - 2*f(i) + f(i-1)) + lambda_f*(f(i)^2 - 1)*f(i) + beta_*s(i)^2*f(i);
 delta(i,2) = (-1/(dx^2))*(s(i+1) - 2*s(i) + s(i-1)) + lambda_s*(s(i)^2 - 1)*s(i) + beta_*f(i)^2*s(i) + ((phase/R)^2)*s(i);
endfor;

delta(1,1) = (-1/(dx^2))*(f(2) - 2*f(1) + f(1)) + lambda_f*(f(1)^2 - 1)*f(1) + beta_*s(1)^2*f(1);
delta(1,2) = (-1/(dx^2))*(s(2) - 2*s(1) + s(1)) + lambda_s*(s(1)^2 - 1)*s(1) + beta_*f(i)^2*s(1) + ((phase/R)^2)*s(1);

delta(n+1,1) = (-1/(dx^2))*(f(n+1) - 2*f(n+1) + f(n)) + lambda_f*(f(n+1)^2 - 1)*f(n+1) + beta_*s(n+1)^2*f(n+1);
delta(n+1,2) = (-1/(dx^2))*(s(n+1) - 2*s(n+1) + s(n)) + lambda_s*(s(n+1)^2 - 1)*s(i+1) + beta_*f(i+1)^2*s(i+1) + ((phase/R)^2)*s(i+1);

grad = dx*[delta(:,1); delta(:,2)];






end
