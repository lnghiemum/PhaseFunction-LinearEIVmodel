function data = datagenSLR(seed, Xtype, n, pY, pW, errortype)
rand('state',seed);
randn('state', seed);

if Xtype==1 % half-normal case
    varX = 1-2/pi;
    X = abs(randn(n,1));
elseif Xtype==2
    % Exponential
    X = exprnd(1,n,1);
    varX = 1;
elseif Xtype==3 %chi-square 3 df
    X = chi2rnd(3,n,1);
    varX = 6;
elseif Xtype==4 %bimodal case
    w = 0.5;
    U1 = rand(n,1);
    X = (U1<w).*normrnd(5,0.6,n,1) + (U1>w).*normrnd(2.5,1,n,1);
    varX = 2.2425; %w = 0.5
end
switch errortype
    case 1 % normal error
        e1 = randn(n,1);
        e2 = randn(n,1);
    case 2 % Laplace error
        e1 =(-log(rand(n,1))+log(rand(n,1)))/sqrt(2);  
        e2 =(-log(rand(n,1))+log(rand(n,1)))/sqrt(2);  
    case 3 % t error with df = 2.5
        e1 = trnd(2.5,n,1)/sqrt(5);
        e2 = trnd(2.5,n,1)/sqrt(5);
end
% beta0 = 0, beta1 = 1 ---> if you change these values, the grid below
% needs to "move" too
beta0 = 1; beta1 = 3;
U = sqrt(pW*varX)*e1;
W = X + U ;
Y = beta0 + beta1*X + sqrt(pY*beta1*varX)*e2;
data = struct('X',X,'Y',Y,'W',W, 'pW',pW); 
end
