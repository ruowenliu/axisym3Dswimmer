function [J, grad] = objectiveDragForceMin(design_vec, THETA_bas, c, compute_grad)
% objective function L_A only for DragForceMin
% Written by Ruowen Liu, 2023

shape = shape3Dadjoint(design_vec);

Cnu = c.multi_cst*(shape.rvol - c.target);

if (min(real(shape.x))>0) && (max(abs(shape.cur))<1e5) % check if x goes negative or curvature too large (very sharp tip)
    J = shape.Jdrag_rByV ...
        - c.lam*Cnu ...
        + c.sig/2*(Cnu^2);
else
    J = Inf;
end

if strcmp(compute_grad,'grad on')
    num_grad = size(THETA_bas,2);
    grad = zeros(num_grad,1);
    for k = 1:num_grad
        shapeDeriv = shape_derivatives(shape, THETA_bas(:,k));
        derivCu = c.multi_cst*shapeDeriv.rvolprime;
        grad(k) = shapeDeriv.Jdrag_rByV_prime - c.lam*derivCu + c.sig*(Cnu)*derivCu;
    end
else
    grad = [];
end

end

function shapeDeriv = shape_derivatives(obj, tht)
% obj is an object (must contain trachat etc.)
% shapeDeriv is a struct
% this is for constrained problem, with struct c

shapeDeriv.tht = tht;
fshat = dotv(obj.trachat, obj.tang);
tht_n = dotv(shapeDeriv.tht,obj.nx);

shapeDeriv.F0prime = -2*pi/obj.mu*(obj.ws'*(fshat.*fshat.*tht_n.*real(obj.x)));

shapeDeriv.Vprime = -2*pi*obj.ws'*(tht_n.*real(obj.x));

Z_deriv = Deriv_panel(imag(obj.x),obj)./obj.sp;
shapeDeriv.Aprime = 2*pi*obj.ws'*( (Z_deriv - (-obj.cur).*real(obj.x)) .* tht_n);

shapeDeriv.rvolprime = 6*sqrt(pi)*( shapeDeriv.Vprime*(obj.area^(-1.5)) -(1.5)*obj.vol*((obj.area)^(-2.5))*shapeDeriv.Aprime );

shapeDeriv.Jdrag_rByV_prime = ((3/4/pi)^(-1/3))/6/pi/obj.mu * ( shapeDeriv.F0prime/((obj.vol)^(1/3)) - obj.F0*shapeDeriv.Vprime/3/((obj.vol)^(4/3)) );

end