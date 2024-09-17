% This file calculates the Bspline representation of the shape
% for a different NL, based on the geometry given.
% Compare the efficiency between NL=21 and NL=99.
% Sep, 2024 @ Ruowen Liu

close all
clear
warningid = 'MATLAB:nearlySingularMatrix'; warning('off',warningid);
addpath('../shape_classes')
addpath('../quadrature_and_kernal')

tic

%% 1. Calculate the geometry based on current NL=21.
file_name = './design_vec.txt';
fID = fopen(file_name, 'r');
design_vec = fscanf(fID, '%f');
fclose(fID);
shape = shape3Dmaxefficiency2(design_vec);
fprintf('Original N_L = 21 \n');
shape.printresults
shapeNL21 = shape; % in the original settings: NL=21

s0.t = shape.t; s0.x = shape.x; 
np = 60;
s0.Z = @(t) interp1(s0.t,real(s0.x),t,'linear','extrap') + 1i*interp1(s0.t,imag(s0.x),t,'linear','extrap');
% After testing, there is no significant differentce among interp1 types.
% s0.Z = @(t) interp1(s0.t,real(s0.x),t,'pchip') + 1i*interp1(s0.t,imag(s0.x),t,'pchip');
% s0.Z = @(t) interp1(s0.t,real(s0.x),t,'spline') + 1i*interp1(s0.t,imag(s0.x),t,'spline');

s0.p = 10; s0.tpan = (0:(pi/np):pi)'; [s0,~] = quadrp(s0, s0.p*np, 'p', 'G'); s0 = arclen(s0);

%% Regenerate the design vector for geometry by B-splines 
NL = 11;
equal_tspace_t_pi = linspace(0,pi,NL+1)';
equal_tspace_x_2pi = [0;real(s0.Z(equal_tspace_t_pi(2:NL)));0;-flip(real(s0.Z(equal_tspace_t_pi(2:NL))));0];
equal_tspace_z_2pi = [imag(s0.Z(0));imag(s0.Z(equal_tspace_t_pi(2:NL)));imag(s0.Z(pi));flip(imag(s0.Z(equal_tspace_t_pi(2:NL))));imag(s0.Z(0))];
MATRIX = zeros(2*NL+5, 2*NL+5);
% dw/dt difference at 0, 2pi
MATRIX(1,1:5) = [-1,-10,0,10,1]/24; MATRIX(1,end-4:end) = -MATRIX(1,1:5);
% d^2w/dt^2 difference at 0, 2pi
MATRIX(2,1:5) = [1,2,-1,2,1]/6; MATRIX(2,end-4:end) = -MATRIX(2,1:5);
% d^3w/dt^3 difference at 0, 2pi
MATRIX(3,1:5) = [-0.5,1,0,-1,0.5]; MATRIX(3,end-4:end) = -MATRIX(3,1:5);
% d^4w/dt^4 difference at 0, 2pi
MATRIX(4,1:5) = [1,-4,6,-4,1]; MATRIX(4,end-4:end) = -MATRIX(4,1:5);
% set up values at exactly equally spaced t locations
for k=5:(2*NL+5)
    % MATRIX is used for solving for xi, note the MATRIX is values of B(t)
    % exactly on equally-spaced grids on [0,L].
    % the four zeros represent equalities of derivatives for periodicity
    MATRIX(k,k-4:k) = [1,26,66,26,1]/120;
end
xi_R_2L= MATRIX\[0;0;0;0;equal_tspace_x_2pi];
xi_Z_2L = MATRIX\[0;0;0;0;equal_tspace_z_2pi];
xi_R = xi_R_2L(1:NL+5); xi_Z = xi_Z_2L(1:NL+5);
design_vec = [xi_R(2:end-1); xi_Z(2:end-1)]; % design vector for shape
shapeNL11 = shape3Dmaxefficiency2_resetNL(NL, design_vec);
fprintf('N_L = %i \n', NL);
shapeNL11.printresults

%% Regenerate the design vector for geometry by B-splines 
NL = 33;
equal_tspace_t_pi = linspace(0,pi,NL+1)';
equal_tspace_x_2pi = [0;real(s0.Z(equal_tspace_t_pi(2:NL)));0;-flip(real(s0.Z(equal_tspace_t_pi(2:NL))));0];
equal_tspace_z_2pi = [imag(s0.Z(0));imag(s0.Z(equal_tspace_t_pi(2:NL)));imag(s0.Z(pi));flip(imag(s0.Z(equal_tspace_t_pi(2:NL))));imag(s0.Z(0))];
MATRIX = zeros(2*NL+5, 2*NL+5);
% dw/dt difference at 0, 2pi
MATRIX(1,1:5) = [-1,-10,0,10,1]/24; MATRIX(1,end-4:end) = -MATRIX(1,1:5);
% d^2w/dt^2 difference at 0, 2pi
MATRIX(2,1:5) = [1,2,-1,2,1]/6; MATRIX(2,end-4:end) = -MATRIX(2,1:5);
% d^3w/dt^3 difference at 0, 2pi
MATRIX(3,1:5) = [-0.5,1,0,-1,0.5]; MATRIX(3,end-4:end) = -MATRIX(3,1:5);
% d^4w/dt^4 difference at 0, 2pi
MATRIX(4,1:5) = [1,-4,6,-4,1]; MATRIX(4,end-4:end) = -MATRIX(4,1:5);
% set up values at exactly equally spaced t locations
for k=5:(2*NL+5)
    % MATRIX is used for solving for xi, note the MATRIX is values of B(t)
    % exactly on equally-spaced grids on [0,L].
    % the four zeros represent equalities of derivatives for periodicity
    MATRIX(k,k-4:k) = [1,26,66,26,1]/120;
end
xi_R_2L= MATRIX\[0;0;0;0;equal_tspace_x_2pi];
xi_Z_2L = MATRIX\[0;0;0;0;equal_tspace_z_2pi];
xi_R = xi_R_2L(1:NL+5); xi_Z = xi_Z_2L(1:NL+5);
design_vec = [xi_R(2:end-1); xi_Z(2:end-1)]; % design vector for shape
shapeNL33 = shape3Dmaxefficiency2_resetNL(NL, design_vec);
fprintf('N_L = %i \n', NL);
shapeNL33.printresults


%% Regenerate the design vector for geometry by B-splines 
NL = 55;
equal_tspace_t_pi = linspace(0,pi,NL+1)';
equal_tspace_x_2pi = [0;real(s0.Z(equal_tspace_t_pi(2:NL)));0;-flip(real(s0.Z(equal_tspace_t_pi(2:NL))));0];
equal_tspace_z_2pi = [imag(s0.Z(0));imag(s0.Z(equal_tspace_t_pi(2:NL)));imag(s0.Z(pi));flip(imag(s0.Z(equal_tspace_t_pi(2:NL))));imag(s0.Z(0))];
MATRIX = zeros(2*NL+5, 2*NL+5);
% dw/dt difference at 0, 2pi
MATRIX(1,1:5) = [-1,-10,0,10,1]/24; MATRIX(1,end-4:end) = -MATRIX(1,1:5);
% d^2w/dt^2 difference at 0, 2pi
MATRIX(2,1:5) = [1,2,-1,2,1]/6; MATRIX(2,end-4:end) = -MATRIX(2,1:5);
% d^3w/dt^3 difference at 0, 2pi
MATRIX(3,1:5) = [-0.5,1,0,-1,0.5]; MATRIX(3,end-4:end) = -MATRIX(3,1:5);
% d^4w/dt^4 difference at 0, 2pi
MATRIX(4,1:5) = [1,-4,6,-4,1]; MATRIX(4,end-4:end) = -MATRIX(4,1:5);
% set up values at exactly equally spaced t locations
for k=5:(2*NL+5)
    % MATRIX is used for solving for xi, note the MATRIX is values of B(t)
    % exactly on equally-spaced grids on [0,L].
    % the four zeros represent equalities of derivatives for periodicity
    MATRIX(k,k-4:k) = [1,26,66,26,1]/120;
end
xi_R_2L= MATRIX\[0;0;0;0;equal_tspace_x_2pi];
xi_Z_2L = MATRIX\[0;0;0;0;equal_tspace_z_2pi];
xi_R = xi_R_2L(1:NL+5); xi_Z = xi_Z_2L(1:NL+5);
design_vec = [xi_R(2:end-1); xi_Z(2:end-1)]; % design vector for shape
shapeNL55 = shape3Dmaxefficiency2_resetNL(NL, design_vec);
fprintf('N_L = %i \n', NL);
shapeNL55.printresults

%% Regenerate the design vector for geometry by B-splines 
NL = 99;
equal_tspace_t_pi = linspace(0,pi,NL+1)';
equal_tspace_x_2pi = [0;real(s0.Z(equal_tspace_t_pi(2:NL)));0;-flip(real(s0.Z(equal_tspace_t_pi(2:NL))));0];
equal_tspace_z_2pi = [imag(s0.Z(0));imag(s0.Z(equal_tspace_t_pi(2:NL)));imag(s0.Z(pi));flip(imag(s0.Z(equal_tspace_t_pi(2:NL))));imag(s0.Z(0))];
MATRIX = zeros(2*NL+5, 2*NL+5);
% dw/dt difference at 0, 2pi
MATRIX(1,1:5) = [-1,-10,0,10,1]/24; MATRIX(1,end-4:end) = -MATRIX(1,1:5);
% d^2w/dt^2 difference at 0, 2pi
MATRIX(2,1:5) = [1,2,-1,2,1]/6; MATRIX(2,end-4:end) = -MATRIX(2,1:5);
% d^3w/dt^3 difference at 0, 2pi
MATRIX(3,1:5) = [-0.5,1,0,-1,0.5]; MATRIX(3,end-4:end) = -MATRIX(3,1:5);
% d^4w/dt^4 difference at 0, 2pi
MATRIX(4,1:5) = [1,-4,6,-4,1]; MATRIX(4,end-4:end) = -MATRIX(4,1:5);
% set up values at exactly equally spaced t locations
for k=5:(2*NL+5)
    % MATRIX is used for solving for xi, note the MATRIX is values of B(t)
    % exactly on equally-spaced grids on [0,L].
    % the four zeros represent equalities of derivatives for periodicity
    MATRIX(k,k-4:k) = [1,26,66,26,1]/120;
end
xi_R_2L= MATRIX\[0;0;0;0;equal_tspace_x_2pi];
xi_Z_2L = MATRIX\[0;0;0;0;equal_tspace_z_2pi];
xi_R = xi_R_2L(1:NL+5); xi_Z = xi_Z_2L(1:NL+5);
design_vec = [xi_R(2:end-1); xi_Z(2:end-1)]; % design vector for shape
shapeNL99 = shape3Dmaxefficiency2_resetNL(NL, design_vec);
fprintf('N_L = %i \n', NL);
shapeNL99.printresults

%%
% shapeNL21.plotblack;
% 
% shapeNL43.plotblack;
% 

%%
toc

d = zeros(4,2);

d(1,1) = abs(shapeNL11.JE - shapeNL99.JE);
d(1,2) = abs(shapeNL11.Jdrag_rByV - shapeNL99.Jdrag_rByV);

d(2,1) = abs(shapeNL21.JE - shapeNL99.JE);
d(2,2) = abs(shapeNL21.Jdrag_rByV - shapeNL99.Jdrag_rByV);

d(3,1) = abs(shapeNL33.JE - shapeNL99.JE);
d(3,2) = abs(shapeNL33.Jdrag_rByV - shapeNL99.Jdrag_rByV);

d(4,1) = abs(shapeNL55.JE - shapeNL99.JE);
d(4,2) = abs(shapeNL55.Jdrag_rByV - shapeNL99.Jdrag_rByV);

fprintf('\n Max difference in E for N_L = 21 and N_L = 99: %e\n', d(2,1));
fprintf('\n Max difference in Jdrag for N_L = 21 and N_L = 99: %e\n', d(2,2));

figure
pA = plot([11,21,33,55], d(:,1),'*--', 'Linewidth', 1.5); hold on; grid on;
pB = plot([11,21,33,55], d(:,2),'o--', 'Linewidth', 1.5); hold on; grid on;
set(gca,'fontsize',17,'TickLabelInterpreter','latex')

ylim([0, 0.021]);
set(gca, 'XTick', [11,21,33,55], ...
    'TickLabelInterpreter','latex');
% set(gca, 'YTick', [1e-4,2e-4,3e-4,4e-4], ...
%     'YTickLabel', {'$1\!\!\times\!\!10^{-4}$','$2\!\!\times\!\!10^{-4}$','$3\!\!\times\!\!10^{-4}$','$4\!\!\times\!\!10^{-4}$'}, ...
%     'TickLabelInterpreter','latex');
xlabel('$N_L$','Interpreter','latex');
legend([pA,pB], ...
    {'$|\Delta E|$ compared to $N_L=99$', ...
    '$|\Delta J_{drag}|$ compared to $N_L=99$'}, ...
    'Interpreter','latex');
