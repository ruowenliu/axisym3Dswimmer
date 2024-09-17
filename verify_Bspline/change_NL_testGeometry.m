% This file calculates the Bspline representation of the shape
% for a different NL, based on the geometry given.
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
shape = shape3Dbasic(design_vec);

s0.t = shape.t; s0.x = shape.x; 
np = 60;
s0.Z = @(t) interp1(s0.t,real(s0.x),t,'linear','extrap') + 1i*interp1(s0.t,imag(s0.x),t,'linear','extrap');
% After testing, there is no significant differentce among interp1 types.
% s0.Z = @(t) interp1(s0.t,real(s0.x),t,'pchip') + 1i*interp1(s0.t,imag(s0.x),t,'pchip');
% s0.Z = @(t) interp1(s0.t,real(s0.x),t,'spline') + 1i*interp1(s0.t,imag(s0.x),t,'spline');

s0.p = 10; s0.tpan = (0:(pi/np):pi)'; [s0,~] = quadrp(s0, s0.p*np, 'p', 'G'); s0 = arclen(s0);

shapeNL21 = shape; % in the original settings: NL=21

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
shapeNL11 = shape3Dbasic_resetNL(NL, design_vec);

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
shapeNL33 = shape3Dbasic_resetNL(NL, design_vec);

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
shapeNL55 = shape3Dbasic_resetNL(NL, design_vec);

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
shapeNL99 = shape3Dbasic_resetNL(NL, design_vec);

%%
linewidth = 1;
markersize = 3;
%%% entire view of R(t)
subplot(2,2,1)
h11 = plot(shapeNL11.t, real(shapeNL11.x),'-', 'LineWidth', linewidth); hold on; grid on;
h21 = plot(shapeNL21.t, real(shapeNL21.x),'-', 'LineWidth', linewidth);
h33 = plot(shapeNL33.t, real(shapeNL33.x),'-', 'LineWidth', linewidth);
h55 = plot(shapeNL55.t, real(shapeNL55.x),'-', 'LineWidth', linewidth);
h99 = plot(shapeNL99.t, real(shapeNL99.x),'-', 'LineWidth', linewidth);
set(gca,'fontsize',20,'TickLabelInterpreter','latex')
xlim([-0.1,pi+0.1])
set(gca, 'XTick', [0, pi/4, pi/2, 3*pi/4, pi], ...
    'TickLabelInterpreter','latex', ...
    'XTickLabel', {'0', '$\pi/4$', '$\pi/2$', '$3\pi/4$', '$\pi$'});
xlabel('$t$','Interpreter','latex');
ylabel('$R$','Interpreter','latex');
legend([h11, h21, h33, h55, h99], ...
    {'$N_L=11$', '$N_L=21$', '$N_L=33$', '$N_L=55$', '$N_L=99$'}, ...
    'Interpreter','latex');

%%% zoom in near pole view of R(t)
subplot(2,2,2)
h11 = plot(shapeNL11.t, real(shapeNL11.x),'.-', 'LineWidth', linewidth, 'MarkerSize', markersize); hold on; grid on;
h21 = plot(shapeNL21.t, real(shapeNL21.x),'s-', 'LineWidth', linewidth, 'MarkerSize', markersize);
h33 = plot(shapeNL33.t, real(shapeNL33.x),'d-', 'LineWidth', linewidth, 'MarkerSize', markersize);
h55 = plot(shapeNL55.t, real(shapeNL55.x),'o-', 'LineWidth', linewidth, 'MarkerSize', markersize);
h99 = plot(shapeNL99.t, real(shapeNL99.x),'*-', 'LineWidth', linewidth, 'MarkerSize', markersize);
set(gca,'fontsize',20,'TickLabelInterpreter','latex')
xlim([2.5, 3.2])
set(gca, 'XTick', [2.6:0.1:3.1, pi], ...
    'TickLabelInterpreter','latex', ...
    'XTickLabel', {'2.6', '2.7', '2.8', '2.9', '3', '3.1', '$\pi$'})
xlabel('$t$','Interpreter','latex');
ylabel('$R$','Interpreter','latex');
legend([h11, h21, h33, h55, h99], ...
    {'$N_L=11$', '$N_L=21$', '$N_L=33$', '$N_L=55$', '$N_L=99$'}, ...
    'Interpreter','latex');

%%% entire view of Z(t)
subplot(2,2,3)
h11 = plot(shapeNL11.t, imag(shapeNL11.x),'-', 'LineWidth', linewidth); hold on; grid on;
h21 = plot(shapeNL21.t, imag(shapeNL21.x),'-', 'LineWidth', linewidth);
h33 = plot(shapeNL33.t, imag(shapeNL33.x),'-', 'LineWidth', linewidth);
h55 = plot(shapeNL55.t, imag(shapeNL55.x),'-', 'LineWidth', linewidth);
h99 = plot(shapeNL99.t, imag(shapeNL99.x),'-', 'LineWidth', linewidth);
set(gca,'fontsize',20,'TickLabelInterpreter','latex')
xlim([-0.1,pi+0.1])
set(gca, 'XTick', [0, pi/4, pi/2, 3*pi/4, pi], ...
    'TickLabelInterpreter','latex', ...
    'XTickLabel', {'0', '$\pi/4$', '$\pi/2$', '$3\pi/4$', '$\pi$'});
xlabel('$t$','Interpreter','latex');
ylabel('$Z$','Interpreter','latex');
legend([h11, h21, h33, h55, h99], ...
    {'$N_L=11$', '$N_L=21$', '$N_L=33$', '$N_L=55$', '$N_L=99$'}, ...
    'Interpreter','latex');

%%% zoom in near pole view of Z(t)
subplot(2,2,4)
h11 = plot(shapeNL11.t, imag(shapeNL11.x),'.-', 'LineWidth', linewidth, 'MarkerSize', markersize); hold on; grid on;
h21 = plot(shapeNL21.t, imag(shapeNL21.x),'.-', 'LineWidth', linewidth, 'MarkerSize', markersize);
h33 = plot(shapeNL33.t, imag(shapeNL33.x),'.-', 'LineWidth', linewidth, 'MarkerSize', markersize);
h55 = plot(shapeNL55.t, imag(shapeNL55.x),'.-', 'LineWidth', linewidth, 'MarkerSize', markersize);
h99 = plot(shapeNL99.t, imag(shapeNL99.x),'.-', 'LineWidth', linewidth, 'MarkerSize', markersize);
set(gca,'fontsize',20,'TickLabelInterpreter','latex')
xlim([2.5, 3.2])
set(gca, 'XTick', [2.6:0.1:3.1, pi], ...
    'TickLabelInterpreter','latex', ...
    'XTickLabel', {'2.6', '2.7', '2.8', '2.9', '3', '3.1', '$\pi$'})
xlabel('$t$','Interpreter','latex');
ylabel('$Z$','Interpreter','latex');
legend([h11, h21, h33, h55, h99], ...
    {'$N_L=11$', '$N_L=21$', '$N_L=33$', '$N_L=55$', '$N_L=99$'}, ...
    'Interpreter','latex');

%%
toc

%%
d = zeros(4,2);

d(1,1) = max(abs(real(shapeNL11.x) - real(shapeNL99.x)));
d(1,2) = max(abs(imag(shapeNL11.x) - imag(shapeNL99.x)));

d(2,1) = max(abs(real(shapeNL21.x) - real(shapeNL99.x)));
fprintf('\n Max difference in R(t) for N_L = 21 and N_L = 99: %e\n', d(2,1));

d(2,2) = max(abs(imag(shapeNL21.x) - imag(shapeNL99.x)));
fprintf('\n Max difference in Z(t) for N_L = 21 and N_L = 99: %e\n', d(2,2));

d(1,1) = max(abs(real(shapeNL33.x) - real(shapeNL99.x)));
d(1,2) = max(abs(imag(shapeNL33.x) - imag(shapeNL99.x)));

d(1,1) = max(abs(real(shapeNL55.x) - real(shapeNL99.x)));
d(1,2) = max(abs(imag(shapeNL55.x) - imag(shapeNL99.x)));

figure
pA = plot([11,21,33,55], d(:,1),'*--', 'Linewidth', 1.5); hold on; grid on;
pB = plot([11,21,33,55], d(:,2),'o--', 'Linewidth', 1.5); hold on; grid on;
set(gca,'fontsize',17,'TickLabelInterpreter','latex')
ylim([0, 0.5e-3]);
set(gca, 'XTick', [11,21,33,55], ...
    'TickLabelInterpreter','latex');
set(gca, 'YTick', [1e-4,2e-4,3e-4,4e-4], ...
    'YTickLabel', {'$1\!\!\times\!\!10^{-4}$','$2\!\!\times\!\!10^{-4}$','$3\!\!\times\!\!10^{-4}$','$4\!\!\times\!\!10^{-4}$'}, ...
    'TickLabelInterpreter','latex');
xlabel('$N_L$','Interpreter','latex');
legend([pA,pB], ...
    {'Maximum of $|\Delta R(t)|$ compared to $N_L=99$', ...
    'Maximum of $|\Delta Z(t)|$ compared to $N_L=99$'}, ...
    'Interpreter','latex');