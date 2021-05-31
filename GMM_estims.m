function parms_out = GMM_estims(W,Y,beta_in)
% Generalized Method of Moments using all moments of the form E(W^iY^j)
% with i+j<=3. This method makes better use of all the information than the
% estimators calculated in MOM.parms.

% Calculates initial values for the other parameters
% For some reason, the measurement error variances are estimated on the
% log-scale. I think I ran into some problems where the numerical routine
% was giving me negative estimated variances, so I built in a safeguard
% against that.
b0_in = beta_in(1); b1_in = beta_in(2);
mu_X = 0.5*mean((Y-b0_in)/b1_in)+0.5*mean(W);
CWY = cov(W,Y);
sig2_X = abs(CWY(1,2)/b1_in);
log_sig2_U = log(abs(CWY(1,1) - sig2_X));
log_sig2_e = log(abs(CWY(2,2) - b1_in.^2.*sig2_X));
mu_3 = mean((W-mu_X).^3);

parms_in = [b0_in,b1_in,mu_X,log(sig2_X),log_sig2_U,log_sig2_e,mu_3];

% The rows of Ind are the powers of W and Y in the GMM approach
Ind = [1,0;0,1;2,0;1,1;0,2;3,0;2,1;1,2;0,3];

Ws = W - mean(W);
Ys = Y - mean(Y);

% Empirical Covariance Matrix for GMM
Sigma = zeros(9,9);
for j = 1:9
    for k = 1:9
        Sigma(j,k) = mean(Ws.^(Ind(j,1)+Ind(k,1)).*Ys.^(Ind(j,2)+Ind(k,2))) - mean(Ws.^Ind(j,1).*Ys.^Ind(j,2))*mean(Ws.^Ind(k,1).*Ys.^Ind(k,2));
    end
end

iSigma = inv(Sigma);

options = optimset('Display','off');
parms_out = fminsearch(@(p) GMM_dist(W,Y,Ind,iSigma,p), parms_in,options);
parms_out(4:6) = exp(parms_out(4:6)); % Exponent returns estimated variances to the postive half-line

end

function D = GMM_dist(W,Y,Ind,iSigma,parms)
% Subroutine that calculates the statistic being minimized using fminsearch
% above.

b0 = parms(1);
b1 = parms(2);
mu_X = parms(3);
sig2_X = exp(parms(4));
sig2_U = exp(parms(5));
sig2_e = exp(parms(6));
mu_3 = parms(7);

n = length(W);
Ws = W - mu_X;
Ys = Y - b0 - b1*mu_X;

mu_vec = [0;0;sig2_X+sig2_U;b1*sig2_X;b1^2.*sig2_X+sig2_e;mu_3;b1*mu_3;b1.^2.*mu_3;b1.^3.*mu_3];

A = zeros(9,1);
for j = 1:9
    A(j,1) = sqrt(1/n)*sum(Ws.^Ind(j,1).*Ys.^Ind(j,2)-mu_vec(j));
end

D = A'*iSigma*A;

end

