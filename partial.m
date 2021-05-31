function result = partial(b,t,W,Y,type)
% This function set up the partial derivative of the phase function
% distance with respect to b0 and b1
if nargin < 5
    type = 1;
end

tt2 = max(t);
dt = t(2)-t(1);
wt = omega(t./tt2,type);

n = length(W);
tvec = repmat(t,n,1);
Yvec = repmat(Y,1,length(t));
lhs_real = nansum(cos(tvec.*Yvec),1);
lhs_im = nansum(sin(tvec.*Yvec),1);

Wvec = repmat(W,1,length(t));
Wb = b(1) + W.* b(2);
Wbvec = repmat(Wb,1,length(t));
rhs_im = nansum(sin(tvec.*Wbvec),1);
rhs_real = nansum(cos(tvec.*Wbvec),1);
G = 2*(lhs_im.*rhs_real-lhs_real.*rhs_im);
G0 = sum(G.*(-lhs_im.*rhs_im.*t - lhs_real.* rhs_real.*t).*(wt).*dt);

A=nansum(tvec.*Wvec.*cos(tvec.*Wbvec),1);
B=nansum(tvec.*Wvec.*sin(tvec.*Wbvec),1);

G1 = nansum(G.*(-lhs_im.*B- lhs_real.*A).*(wt).*dt);
result = [G0 G1];
end