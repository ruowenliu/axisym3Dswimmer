function intval = integral_to_quadrp(s, integ)
% integral from 0 to each quadrature point

intval = zeros(size(s.t));

for k = 1:length(s.t)
tEnd = s.t(k); idx = (s.thi<tEnd); 
tlo = s.tlo(sum(idx)+1); thi = tEnd; [x0, w0, D0] = gauss(s.p);
t = tlo+(1+x0)/2*(thi-tlo); w = w0*(thi-tlo)/2; 
x = interpolate_by_panel(s, s.x, t); % x = s.Z(t); 
xp = D0*x*2/(thi-tlo); ws = w.*abs(xp);
F = sum(integ(1:s.p*sum(idx)).*s.ws(1:s.p*sum(idx))); 
L = interpmat_1d(-1+(1+x0)*(thi-tlo)/(s.thi(sum(idx)+1)-s.tlo(sum(idx)+1)),x0); fL = L*integ(s.p*sum(idx)+(1:s.p));
F = F+sum(fL.*ws);
intval(k) = F;
end

end

function L = interpmat_1d(t,s)
p = numel(s); q = numel(t); s = s(:); t = t(:);       % all col vecs
n = p; % set the polynomial order we go up to (todo: check why bad if not p)
V = ones(p,n); for j=2:n, V(:,j) = V(:,j-1).*s; end   % polyval matrix on nodes
R = ones(q,n); for j=2:n, R(:,j) = R(:,j-1).*t; end   % polyval matrix on targs
L = (V'\R')'; % backwards-stable way to do it (Helsing) See corners/interpdemo.m
end