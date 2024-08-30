function shapeDeriv = shape_derivatives(obj, tht)
% obj is an object (must contain trachat etc.)
% shapeDeriv is a struct

shapeDeriv.tht = tht;
fshat = dotv(obj.trachat, obj.tang);
tht_n = dotv(shapeDeriv.tht,obj.nx);

shapeDeriv.F0prime = -2*pi/obj.mu*(obj.ws'*(fshat.*fshat.*tht_n.*real(obj.x)));

shapeDeriv.volprime = -2*pi*obj.ws'*(tht_n.*real(obj.x));

Z_deriv = Deriv_panel(imag(obj.x),obj)./obj.sp;
shapeDeriv.areaprime = 2*pi*obj.ws'*( (Z_deriv - (-obj.cur).*real(obj.x)) .* tht_n);

shapeDeriv.rvolprime = 6*sqrt(pi)*( shapeDeriv.volprime*(obj.area^(-1.5)) -(1.5)*obj.vol*((obj.area)^(-2.5))*shapeDeriv.areaprime );

shapeDeriv.ellprime = - obj.ws'*((- obj.cur).*tht_n);

shapeDeriv.Uprime = ShapeSens_U(obj, shapeDeriv.tht);

shapeDeriv.JWprime = ShapeSens_Jw(obj, shapeDeriv.tht);

shapeDeriv.JDprime = shapeDeriv.F0prime * obj.U * obj.U + obj.F0 * 2 * obj.U * shapeDeriv.Uprime;

shapeDeriv.JEprime = (shapeDeriv.JDprime - obj.JE*shapeDeriv.JWprime)/obj.JW;

%
shapeDeriv.Jdrag_rByVprime = ((3/4/pi)^(-1/3))/6/pi/obj.mu * ( shapeDeriv.F0prime/((obj.vol)^(1/3)) - obj.F0*shapeDeriv.volprime/3/((obj.vol)^(4/3)) );

end

function sensU = ShapeSens_U(s, tht)
% this function computes the shape sensitivity of U
% given shape transformations velocity \theta (tht)

kap = - s.cur;
tht_s = dotv(tht,s.tang);
tht_n = dotv(tht,s.nx);
fs    = dotv(s.trac, s.tang);
fshat = dotv(s.trachat, s.tang);

R = real(s.x);

Dus = Deriv_panel(s.uslip,s)./s.sp;

DRusthtn = Deriv_panel((R.*s.uslip.*tht_n),s)./s.sp;

F0 = 2*pi*s.ws'*(imag(s.trachat).*real(s.x));

integ = kap.*tht_n;
s_star = tht_s - integral_to_quadrp(s, integ);

el = sum(s.ws);
el_star = - s.ws'*(kap.*tht_n);

term1 = R.*fshat.*(s_star-tht_s-s.arclen/el*el_star).*Dus;

term2 = -s.preshat.*DRusthtn;

term3a = kap.*R.*s.uslip.*fshat;

term3b = -1/s.mu*R.*fs.*fshat;

term3 = (term3a+term3b).*tht_n;

sensU = -2*pi/F0*(s.ws'*(term1+term2+term3));

end % sensU

function sensJw = ShapeSens_Jw(s, tht)
% this function computes the shape sensitivity of Jw
% given shape transformations velocity \theta (tht)

kap = - s.cur;
tht_s = dotv(tht,s.tang);
tht_n = dotv(tht,s.nx);
fs    = dotv(s.trac, s.tang);
fn    = dotv(s.trac, s.nx);

Dtht_n = Deriv_panel(tht_n,s)./s.sp;

R = real(s.x);
DR = Deriv_panel(R,s)./s.sp;

Dus = Deriv_panel(s.uslip,s)./s.sp;

DRus = Deriv_panel((R.*s.uslip),s)./s.sp;

integ = kap.*tht_n;
s_star = tht_s - integral_to_quadrp(s, integ);

el = sum(s.ws);
el_star = - s.ws'*(kap.*tht_n);

term1 = R.*fs.*(s_star-tht_s-s.arclen/el*el_star).*Dus;

term2 = s.uslip.*fn.*Dtht_n.*R;

term3a = -2*s.mu*Dus.*DR.*s.uslip;

term3b = kap.*R.*s.uslip.*fs;

term3c = -.5/s.mu*R.*(fs.^2);

term3d = -s.pres.*DRus;

term3 = (term3a+term3b+term3c+term3d).*tht_n;

sensJw = 4*pi*(s.ws'*(term1+term2+term3));

end % sensJW

