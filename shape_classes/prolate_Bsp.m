function [xi_R, xi_Z] = prolate_Bsp(np, axis, NL, displayerrorcheck)
% (c) 2023 Ruowen Liu
% regenerate shape (prolate spheroid) using B-splines
% It uses a full rotation (2pi) to generate initial xi_R, xi_Z

% % prolate spheroid
s0.Z = @(t) axis.a*sin(t)+1i*axis.b*cos(t); 
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
