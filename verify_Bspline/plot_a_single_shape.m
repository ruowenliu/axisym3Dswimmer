% This file tests different settings of Bspline to represent the geometry
% Sep, 2024 @ Ruowen Liu

% IMPORTANT: Must manually open ../shape_classes/shape3Dparam.m and set NL
% {mustBeNumeric} = 21 ot a different number
% Pay

close all
clear
warningid = 'MATLAB:nearlySingularMatrix'; warning('off',warningid);
addpath('../shape_classes')
addpath('../quadrature_and_kernal')

% Read design vector
file_name = './design_vec.txt';
fID = fopen(file_name, 'r');
design_vec = fscanf(fID, '%f');
fclose(fID);
shape = shape3Dmaxefficiency2(design_vec);

shape.U

shape.printresults

%%
shape.plotuslip

saveas(gcf, 'iter17_nu_070_uslip', 'epsc');

%%
shape.plotblue

saveas(gcf, 'iter17_nu_070_shape', 'epsc');

%%

figure
plot(shape.t, real(shape.x), 'LineWidth', 1.5);
xlabel('$t$','Interpreter','latex');
ylabel('$R$','Interpreter','latex');
set(gca, 'fontsize', 20);
xlim([-0.1,pi+0.1])
set(gca, 'XTick', [0, pi/4, pi/2, 3*pi/4, pi], ...
    'TickLabelInterpreter','latex', ...
    'XTickLabel', {'0', '$\pi/4$', '$\pi/2$', '$3\pi/4$', '$\pi$'})
ylim([-0.1,0.5])
grid on;
drawnow;

saveas(gcf, 'iter17_nu_070_R', 'epsc');

%%

figure
plot(shape.t, imag(shape.x), 'LineWidth', 1.5);
xlabel('$t$','Interpreter','latex');
ylabel('$Z$','Interpreter','latex');
set(gca, 'fontsize', 20);
xlim([-0.1,pi+0.1])
set(gca, 'XTick', [0, pi/4, pi/2, 3*pi/4, pi], ...
    'TickLabelInterpreter','latex', ...
    'XTickLabel', {'0', '$\pi/4$', '$\pi/2$', '$3\pi/4$', '$\pi$'})
ylim([-2.1,2.1])
grid on;
drawnow;

saveas(gcf, 'iter17_nu_070_Z', 'epsc');


