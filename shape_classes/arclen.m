function s = arclen(s)
% arc length up to each quadrature point
% last modified on May 14, 2022 by Ruowen Liu

s.arclen = zeros(size(s.t));

for k = 1:length(s.t)
tEnd = s.t(k); idx = (s.thi<tEnd); arcLen = sum(s.ws(1:s.p*sum(idx))); % entire panel \subset [0,t]; 
tlo = s.tlo(sum(idx)+1); thi = tEnd; [x0, w0, D0] = gauss(s.p); 
t = tlo+(1+x0)/2*(thi-tlo); w = w0*(thi-tlo)/2; 
x = s.Z(t); xp = D0*x*2/(thi-tlo);
ws = w.*abs(xp);
s.arclen(k) = arcLen+sum(ws);
end

s.ell = sum(s.ws);
s.vol = body_volume(s);
s.area = body_area(s);
s.rvol = 6*sqrt(pi)*s.vol/(s.area^(3/2));
end