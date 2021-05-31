function est = estimatorSLR(Y,W)
% the function computes estimators for simple linear EIV model
% Y : n*1 vector of response
% W : n*1 vector of contaminated covariates
% type: weight function used in computing the phase function distance
n = length(Y);
%% Naive estimator
Wmat = [ones(length(W),1) W];
beta_naive = Wmat\Y;

%% GMM estimator;
beta_in = MOM_parms(W,Y);
beta_GMM = GMM_estims(W,Y,beta_in);
beta_GMM_est = beta_GMM(1:2);

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
t = linspace(0.0001*t_star,t_star,1001);

options = optimoptions('fsolve','Display','off');
f1= @(b)partial(b,t,W,Y);
b_phase1= fsolve(f1,beta_GMM_est,options);

est = struct('bnaive',beta_naive','bGMM',beta_GMM_est,'bPhase',b_phase1, 'tGrid',t);
end



