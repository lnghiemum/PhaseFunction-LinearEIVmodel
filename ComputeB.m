function B = ComputeB(W , Y , b , t, type)

tt2 = max(t);
wt = omega(t./tt2, type);
dt = t(2)-t(1);

n = length(W);
tvec = repmat(t,n,1);
Yvec = repmat(Y,1,length(t));
Wvec = repmat(W,1,length(t));
Wb = b(1) + b(2).*W;
Wbvec = repmat(Wb,1,length(t));

% Caculate sine and cosine components at t-values for numerical integration
cos_Y = sum(cos(tvec.*Yvec),1);
sin_Y = sum(sin(tvec.*Yvec),1);
sin_Wb = sum(sin(tvec.*Wbvec),1);
cos_Wb = sum(cos(tvec.*Wbvec),1);
Wsin_Wb = sum(Wvec.*sin(tvec.*Wbvec),1);
Wcos_Wb = sum(Wvec.*cos(tvec.*Wbvec),1);
WWsin_Wb = sum(Wvec.^2 .*sin(tvec.*Wbvec),1);
WWcos_Wb = sum(Wvec.^2 .*cos(tvec.*Wbvec),1);

% Components of B matrix
I_11 = dt*sum(2*t.^2.*wt.*(cos_Y.*cos_Wb+sin_Y.*sin_Wb).^2 - ...
    2*t.^2.*wt.*(sin_Y.*cos_Wb-cos_Y.*sin_Wb).^2);

I_12 = dt*sum(2*t.^2.*wt.*(cos_Y.*cos_Wb+sin_Y.*sin_Wb).*(cos_Y.*Wcos_Wb+sin_Y.*Wsin_Wb) - ...
    2*t.^2.*wt.*(sin_Y.*cos_Wb-cos_Y.*sin_Wb).*(sin_Y.*Wcos_Wb-cos_Y.*Wsin_Wb));

I_22 = dt*sum(2*t.^2.*wt.*(cos_Y.*Wcos_Wb+sin_Y.*Wsin_Wb).^2 - ...
    (2*t.^2.*wt.*(cos_Y.*sin_Wb - sin_Y.*cos_Wb) .* ...
    (cos_Y.*WWsin_Wb - sin_Y.*WWcos_Wb)));

B = [I_11,I_12;I_12,I_22];