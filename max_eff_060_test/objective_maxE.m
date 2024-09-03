function [J, grad] = objective_maxE(new_design_vec,THETA_bas,c,compute_grad,fixed_RN,fixed_RS,fixed_ZN,fixed_ZS)
% objective function L_A
% design_vec reduced dim

design_vec = fixpolesreverse(new_design_vec,fixed_RN,fixed_RS,fixed_ZN,fixed_ZS);

shape = shape3Dmaxefficiency2(design_vec);

Cnu = c.multi_cst*(shape.rvol - c.target);

if (min(real(shape.x))>0) && (max(abs(shape.cur))<1e5) && (shape.JE>0) % check if x goes negative or curvature too large (very sharp tip) or Max-eff is positive
    J = - shape.JE ...
        - c.lam*(Cnu) ...
        + c.sig/2*(Cnu^2);
else
    J = Inf;
end

fixed_dim = length(fixed_RN);
new_half_dim = length(new_design_vec)/2;

if strcmp(compute_grad,'grad on')
    grad = zeros(length(new_design_vec),1);
    for k = 1:new_half_dim
        shapeDer = shape_derivatives(shape, THETA_bas(:,k+fixed_dim));
        % grad for R portion
        grad(k) = - shapeDer.JEprime ...
            - c.lam*c.multi_cst*shapeDer.rvolprime + c.sig*(Cnu)*c.multi_cst*shapeDer.rvolprime;
        % grad for Z portion
        shapeDer = shape_derivatives(shape, THETA_bas(:,k+new_half_dim+3*fixed_dim));
        grad(k+new_half_dim) = - shapeDer.JEprime ...
            - c.lam*c.multi_cst*shapeDer.rvolprime + c.sig*(Cnu)*c.multi_cst*shapeDer.rvolprime;
    end % for 1
else
    grad = [];
end

end

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

shapeDeriv.Uprime = ShapeSens_U(obj, shapeDeriv.tht);

shapeDeriv.JWprime = ShapeSens_Jw(obj, shapeDeriv.tht);

shapeDeriv.JDprime = shapeDeriv.F0prime * obj.U * obj.U + obj.F0 * 2 * obj.U * shapeDeriv.Uprime;

shapeDeriv.JEprime = (shapeDeriv.JDprime - obj.JE*shapeDeriv.JWprime)/obj.JW;

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

