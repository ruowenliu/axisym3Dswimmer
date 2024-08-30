function design_vec = get_initial_parameterized_shape_gallery(shapetype, nu)
% Copyright @ Ruowen Liu, June 2023

if strcmp(shapetype,"prolate")
    design_vec = get_initial_prolate(nu, "no noise");
else
    shape = shape3Dparam; % Change setting in "shape3Dparam": np and NL
    [xi_R_2L, xi_Z_2L] = Bsp_parameterized_shape(shapetype, shape.np, shape.NL);
    xi_R = xi_R_2L(1:shape.NL+5); xi_Z = xi_Z_2L(1:shape.NL+5);
    design_vec = [xi_R(2:end-1); xi_Z(2:end-1)]; % design vector for shape
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xi_R, xi_Z] = Bsp_parameterized_shape(shapetype, np, NL, displayerrorcheck)

switch shapetype
    case "symmetric-peanut-large"
        ca = 0.25; cb = 2-2*ca; cc = 1.6; cd = 1.8; ce = 0.15;
        s0.Z = @(t) 1.2 * ( ((cos(cb*t+ca*pi)).^5-ce*(cos(cb*t+ca*pi)).^2-(cos(cb*t+ca*pi))-(cos(ca*pi)).^5+ce*(cos(ca*pi)).^2+(cos(ca*pi)))/cc + 1i*cd*cos(t) );
    case "symmetric-peanut-small"
        ca = 0.25; cb = 2-2*ca; cc = 1.6; cd = 1.8; ce = 0.15;
        s0.Z = @(t) 0.8 * ( ((cos(cb*t+ca*pi)).^5-ce*(cos(cb*t+ca*pi)).^2-(cos(cb*t+ca*pi))-(cos(ca*pi)).^5+ce*(cos(ca*pi)).^2+(cos(ca*pi)))/cc + 1i*cd*cos(t) );
    case "symmetric-peanut"
        ca = 0.25; cb = 2-2*ca; cc = 1.6; cd = 1.8; ce = 0.15;
        s0.Z = @(t) ((cos(cb*t+ca*pi)).^5-ce*(cos(cb*t+ca*pi)).^2-(cos(cb*t+ca*pi))-(cos(ca*pi)).^5+ce*(cos(ca*pi)).^2+(cos(ca*pi)))/cc + 1i*cd*cos(t);
    case "asymmetric-vase"
        ca = 0.25; cb = 2-2*ca; cc = 1.6; cd = 1.6; ce = 0.15;
        s0.Z = @(t) ((cos(cb*t+ca*pi)).^5-ce*(cos(cb*t+ca*pi)).^2-(cos(cb*t+ca*pi))-(cos(ca*pi)).^5+ce*(cos(ca*pi)).^2+(cos(ca*pi)))/cc + 1i*cd*cos(1/pi*t.^2);
    case "three-bumps" % This is a problematic shape. May need to increase discretization
        xfunc = @(ca,x) 0.001*(-(5*x*ca-1).*(x*ca-1).*(x*ca-pi).*(x*ca-4).*(x*ca-5).*(x*ca-6.8)+408);
        afun = @(ca) xfunc(ca,pi);
        ca_val = fzero(afun,2.203);
        s0.Z = @(t) xfunc(ca_val,t) + 1i * 1.65 * cos(t);
    case "three-bumps-2"
        xfunc = @(ca,x) 0.001*(-(3.5*x*ca-1).*(x*ca-1).*(x*ca-pi).*(x*ca-4).*(x*ca-5).*(x*ca-6.8)+408);
        afun = @(ca) xfunc(ca,pi);
        ca_val = fzero(afun,2.203);
        s0.Z = @(t) (xfunc(ca_val,t)).*((-t).*(t-pi)) + 1i * 1.7 * cos(t);
    case "long-lump"
        s0.Z = @(t) -0.1*(t.*(t-pi).*(0.6*t.^3-0.9*t.^2-1.2*t+4))+ 1i * 1.8 * cos(1/pi*t.^2);
    case "peanut-small-head"
        cf = 0.7315;
        s = @(t) 1-cf*cos(pi*(cos(t)-.2)/2);
        s0.Z = @(t) (sqrt((1-cos(t).^2).*s(t).^2) + 1i*cos(t)) * 1.5;
    case "peanut-round"
        cf = 0.7;
        s = @(t) 1-cf*cos(0.5*pi*cos(t));
        s0.Z = @(t) (sqrt((1-cos(t).^2)).*s(t) + 1i*cos(t)) * 1.5;
    case "peanut-bigger"
        cst1 = 1.5;
        cst2 = cst1*(cos(0.3*pi))^2+cos(0.3*pi);
        s0.Z = @(t) -cst1*(cos(1.4*t+0.3*pi).^2)-cos(1.4*t+0.3*pi)+cst2 + 1i * 1.6 * cos(t);
end

s0.p = 10; s0.tpan = (0:(pi/np):pi)'; [s0,~] = quadrp(s0, s0.p*np, 'p', 'G'); s0 = arclen(s0);

% % regenerate by B-splines
equal_tspace_t_pi = linspace(0,pi,NL+1)';
equal_tspace_x_2pi = [0;real(s0.Z(equal_tspace_t_pi(2:NL)));0;-flip(real(s0.Z(equal_tspace_t_pi(2:NL))));0];
equal_tspace_z_2pi = [imag(s0.Z(0));imag(s0.Z(equal_tspace_t_pi(2:NL)));imag(s0.Z(pi));flip(imag(s0.Z(equal_tspace_t_pi(2:NL))));imag(s0.Z(0))];
MATRIX = zeros(2*NL+5, 2*NL+5);
% % dw/dt difference at 0, 2pi
MATRIX(1,1:5) = [-1,-10,0,10,1]/24; MATRIX(1,end-4:end) = -MATRIX(1,1:5);
% % d^2w/dt^2 difference at 0, 2pi
MATRIX(2,1:5) = [1,2,-1,2,1]/6; MATRIX(2,end-4:end) = -MATRIX(2,1:5);
% % d^3w/dt^3 difference at 0, 2pi
MATRIX(3,1:5) = [-0.5,1,0,-1,0.5]; MATRIX(3,end-4:end) = -MATRIX(3,1:5);
% % d^4w/dt^4 difference at 0, 2pi
MATRIX(4,1:5) = [1,-4,6,-4,1]; MATRIX(4,end-4:end) = -MATRIX(4,1:5);
% % set up values at exactly equally spaced t locations
for k=5:(2*NL+5), MATRIX(k,k-4:k) = [1,26,66,26,1]/120; end
% % MATRIX is used for solving for xi, note the MATRIX is values of B(t)
% % exactly on equally-spaced grids on [0,L].
% % the four zeros represent equalities of derivatives for periodicity
xi_R = MATRIX\[0;0;0;0;equal_tspace_x_2pi];
xi_Z = MATRIX\[0;0;0;0;equal_tspace_z_2pi];

% % verify the error of s0 and s
if nargin > 3
    disp(displayerrorcheck);
    L = pi;
    fun = @(t) B0_original( (t*(2*NL)/(2*L)) - ((-5):((2*NL)-1)) ) * ( xi_R + 1i * xi_Z );
    t1 = linspace(0,pi,500)'; disp(max(abs(s0.Z(t1)-fun(t1))));
end
end