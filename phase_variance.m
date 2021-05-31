function res = phase_variance(Y,W, tGrid, betaHat, B, full_bootstrap, type)
% The function computes variance of phase function using two bootstrap
% procedures
% B = number of bootstrap samples 
% tGrid: the grid of t used in computing betaHat
% if you want full bootstrap algorithm (more computationally expensive), 
% you can set full_bootstrap = 1 
if nargin==5
    type = 1;
    full_bootstrap = 0;   
elseif nargin==6
    type =1;
end
    
nObs = length(Y);
%% Modified Bootstrap algorithm
G_out=zeros(B,2);
for b = 1:B
index = ceil(nObs*rand(nObs,1));
Wb = W(index); Yb = Y(index);
n=size(Wb,1);
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
t = linspace(0.0001*t_star,t_star,1001);
G_out(b,:) = partial(betaHat,t,Wb,Yb,type);
end

B = ComputeB(W , Y , betaHat , tGrid, type);
A = cov(G_out);

plug_in_var = inv(B)*A*inv(B);

res = struct('plug_in_var',plug_in_var);
%% Full Bootstrap
if full_bootstrap == 1
    boot_var = bootstrap(W,Y, B);
    res = struct('plug_in_var',plug_in_var, 'full_bootstrap',boot_var);
end 

end


