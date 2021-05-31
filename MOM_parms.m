function [beta_out3,beta_out4,beta_out5] = MOM_parms(W,Y)
% This calculates (1)beta_out3, an estimator based on the skewness of W and
% Y, (2) beta_out4, an estimator based on Cov(W,Y^2) and Cov(W^2,Y) and (3)
% an estimator that averages out the other two estimators

% None of these are of primary interest, but are used to find initial
% values for the Generalized Method of Moments algorithm.
b1_out3 = nthroot(mean((Y-mean(Y)).^3),3)/nthroot(mean((W-mean(W)).^3),3);
b0_out3 = mean(Y)-b1_out3*mean(W);
beta_out3 = [b0_out3,b1_out3];
b1_out4 = mean((W-mean(W)).*(Y-mean(Y)).^2)/mean((W-mean(W)).^2.*(Y-mean(Y)));
b0_out4 = mean(Y)-b1_out4*mean(W);
beta_out4 = [b0_out4,b1_out4];

b0 = (b0_out3+b0_out4)/2;
b1 = (b1_out3+b1_out4)/2;
beta_out5 = [b0,b1];

end