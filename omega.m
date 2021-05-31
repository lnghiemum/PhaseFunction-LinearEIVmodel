function wt = omega(t,type, delta)
% This function computes the  function K(t) used in the formula for the
% phase function estimator
    if nargin==3
        f = @(t) exp(-delta/2*t.^2);
    elseif type==1
        f = @(t) (1-t).^2;
    elseif type==2
        f = @(t) (1-abs(t));
    elseif type==3
        f = @(t) (1-t.^(1/2));
    elseif type==4
        f = @(t) (1-t.^2);    
    end
    wt = f(t);
end