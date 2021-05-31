function [X,Z] = datagenCopula(seed,Xtype, rho, n)
% The function uses normal copula model to generate X and Z of any distribution
% with correlation rho
rng(seed);
C = [1,rho;rho,1];
Z0 = chol(C)'*randn(2,n);
Z0 = Z0';
U0 = normcdf(Z0);
options = optimoptions('fsolve','Display','off');

if Xtype==1 % half-normal case
    for k = 1:n
        X(k,1) = fsolve(@(x) 2*normcdf(sqrt(2)*x)-1-U0(k,1),1,options);
        Z(k,1) = fsolve(@(x) 2*normcdf(sqrt(2)*x)-1-U0(k,2),1,options);
    end
elseif Xtype==4 %bimodal case
    for k = 1:n
        X(k,1) = fsolve(@(x) 0.5*normcdf((x-5)/0.6)+0.5*normcdf((x-2.5)/1)-U0(k,1),4,options);
        Z(k,1) = fsolve(@(x) 0.5*normcdf((x-5)/0.6)+0.5*normcdf((x-2.5)/1)-U0(k,2),4,options);
    end
end