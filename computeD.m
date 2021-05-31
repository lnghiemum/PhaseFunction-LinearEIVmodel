function D = computeD(b,t,W,Y,type)
% The function computes the value of the D statistic
% b is the value of coefficients, writing as a COLUMN vector
% t is the vector of grid values of t for approximating the integral
% W is a matrix of covariates
% Y is the response
% type is the type of weighting function

tt2 = max(t);
dt = t(2)-t(1);
wt = omega(t./tt2,type);

n = length(Y);
tvec = repmat(t,n,1);
Yvec = repmat(Y,1,length(t));
%phi_Y = mean(exp(tvec.*Yvec*1i),1);
lhs_real = sum(cos(tvec.*Yvec),1);
lhs_im = sum(sin(tvec.*Yvec),1);

Wb = W*b;
Wbvec = repmat(Wb,1,length(t));
rhs_im = sum(sin(tvec.*Wbvec),1);
rhs_real = sum(cos(tvec.*Wbvec),1);

D = sum((lhs_real.*rhs_im-lhs_im.*rhs_real).^2 .* wt*dt);
D = log(D);

end

