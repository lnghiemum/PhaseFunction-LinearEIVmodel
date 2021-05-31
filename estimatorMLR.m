function est = estimatorMLR(Y,W)
% compute the phase function estimator for MLR EIV model. 
% Y: outcome vector
% W: matrix with both error-prone and error-free covariates

% both measured with error and without error
n = length(Y);
p = size(W,2);
W = W-repmat(mean(W),n,1);  % center W
Y = Y-mean(Y); % center Y
%% Naive estimator
beta_naive = W\Y;

%% Phase Function Regression Estimators
% Find the t_star value for the phase function interval
ref_val = n^(-0.25);
k = 0;
t = linspace(k,k+1,50);
tvec = repmat(t,n,1);
Yvec = repmat(Y,1,50);
phi_t = mean(exp(1i*tvec.*Yvec));
mod_phi = sqrt(real(phi_t.*conj(phi_t)));
while min(mod_phi)>ref_val
    k = k+1;
    t = linspace(k,k+1,50);
    tvec = repmat(t,n,1);
    phi_t = mean(exp(1i*tvec.*Yvec));
    mod_phi = sqrt(real(phi_t.*conj(phi_t)));
end
t_star = interp1(mod_phi,t,ref_val);
t = linspace(0.001*t_star,t_star,1001);
b = beta_naive;

f = @(b) computeD(b,t,W,Y,1);
b_phase = fminsearch(f,b);
est = struct('b_naive',beta_naive,'b_phase',b_phase);
end



