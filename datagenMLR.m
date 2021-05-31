function data = datagenMLR(seed, Xtype, n, pY, pW, errortype, rho)
% This function generates data used in the simulation of the MLR EIV model
[X,Z]=datagenCopula(seed, Xtype, rho, n);
if Xtype==1 
    varX =1-2/pi;
elseif Xtype==4 
    varX = 2.2425;
end    
switch errortype
    case 1
        e1 = randn(n,1);
        e2 = randn(n,1);
    case 2
        e1 =(-log(rand(n,1))+log(rand(n,1)))/sqrt(2);  
        e2 =(-log(rand(n,1))+log(rand(n,1)))/sqrt(2);  
    case 3
        e1 = trnd(2.5,n,1)/sqrt(5);
        e2 = trnd(2.5,n,1)/sqrt(5);
end
% beta0 = 0, beta1 = 1 ---> if you change these values, the grid below
% needs to "move" too
beta0 = 0; beta1 = 3; beta2=2;
U = sqrt(pW*varX)*e1;
W = X + U ;
Y = beta0+ beta1*X + beta2*Z+ sqrt(pY*beta1*varX)*e2;
data = struct('X',X,'Y',Y,'W',W, 'pW',pW, 'Z',Z); 
end
