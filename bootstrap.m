function res = bootstrap(W,Y,B)
% The function computes the variance for the phase function estimator using
% the full bootstrap procedure for the EIV simple linear model
% W: contaminated covariate; Y: response; B: number of bootstrap samples

n = length(Y);
for b = 1:B
    index = ceil(n*rand(n,1));
    Wb = W(index); Yb = Y(index);
    beta_in = MOM_parms(Wb,Yb);
    beta_GMM = GMM_estims(W,Y,beta_in);
    beta_GMM_est = beta_GMM(1:2);
    ref_val = n^(-0.25);
    k = 0;
    t = linspace(k,k+1,50);
    tvec = repmat(t,n,1);
    Yvec = repmat(Yb,1,50);
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
    t = linspace(0,t_star,1001);
    options = optimoptions('fsolve','Display','off');
    f1= @(b)partial(b,t,Wb,Yb,1);
    [b_phase1,~,out] = fsolve(f1,beta_GMM_est,options);
    phase_boot(b,:) = b_phase1;
end
res =[cov(phase_boot)].*(B-1)/B;
end

