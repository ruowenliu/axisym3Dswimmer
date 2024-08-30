function fp = Deriv_panel(fval,s)
% This function computes the first derivative w.r.t. t using panels
% Arguments:
%
%   fval -- a function defined on the wall
%
%   s -- structure of the wall
%
% Return:
%
%   fp  --   the first derivative
%

p = s.p;
pt = s.thi - s.tlo;
np = s.np;

[~, ~, D] = gauss(p);

fp = zeros(length(fval),1);

for i=1:np
    ii = (i-1)*p+(1:p); % indices of this panel
    fp(ii) = D*fval(ii)*2/pt(i);
end

end