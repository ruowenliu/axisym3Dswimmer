%%% preset
clear; close all; format long; format compact;
warningid = 'MATLAB:nearlySingularMatrix';
warning('off',warningid)
addpath('../shape_classes')
addpath('../quadrature_and_kernal')

%%% read the shape with minimum drag force
nu = 0.70; c.sig = 1.2 * 1e4;
dir_name = '../min_drag_force_various_nu/';
file_name = [dir_name 'min_drag_main_nu_' num2str(100*nu, '%.3i') '.mat'];
load(file_name, 'design_vec');
design_vec_0 = design_vec;
shape_initial = shape3Dmaxefficiency2(design_vec_0);
shape_initial.plotgreen;

tStart = tic;

c.target = nu;
c.lam = 0;
fprintf('* Reset c.sig %g \n',c.sig);
c.multi_cst = 1;

fprintf('--- Initial ---\n');
shape_initial = shape3Dmaxefficiency2(design_vec_0);
shape_initial.printresults;
shape_initial.plotgreen;
% saveas(gcf, ['./' resultdir '/shape_'  num2str(nu*100,'%2i') '_initial_dragmin'],'pdf')
% saveas(gcf, ['./' resultdir '/shape_'  num2str(nu*100,'%2i') '_initial_dragmin'],'epsc')
THETA_bas = matrix_theta_basis(shape_initial.t, shape_initial.NL, shape_initial.L);

options = optimoptions(@fminunc, ...
    'Display','iter', ...
    'Algorithm','quasi-newton', ...
    'HessianApproximation','lbfgs', ...
    'SpecifyObjectiveGradient',true, ...
    'UseParallel',true ...    
    ); 
%    'OutputFcn',@outfunMaxEffi, ...
% 'StepTolerance', 1e-3 ...

L_A = @(vec) objectiveMaxEffi(vec, THETA_bas, c, 'grad on');
[design_vec_final,~,~,output] = fminunc(L_A,design_vec_0,options);

shape_final = shape3Dmaxefficiency2(design_vec_final);
shape_final.printresults;
shape_final.plotorange;

% save(['./' resultdir '/data_' num2str(nu*100,'%2i') '_final.mat']);

% saveas(gcf, ['./' resultdir '/shape_'  num2str(nu*100,'%2i') '_final_maxeffi'],'pdf')
% saveas(gcf, ['./' resultdir '/shape_'  num2str(nu*100,'%2i') '_final_maxeffi'],'epsc')

tEnd = toc(tStart);
fprintf('Total Elapsed Time: %i hours, %i minutes, %i seconds. \n', ...
    round(tEnd/3600), round(mod(tEnd,3600)/60), ceil(mod(tEnd,60)));

diary off

