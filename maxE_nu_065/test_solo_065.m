% This file produces the shape optimization process to maximize
% swimming efficiency
% constraned with reduced volume nu
% Written by Ruowen Liu, July 2023
% Last modified and tested in September 2024

%%% preset
clear; close all; format long; format compact;
warningid = 'MATLAB:nearlySingularMatrix';
warning('off',warningid)
addpath('../shape_classes')
addpath('../quadrature_and_kernal')
addpath('../optimization_functions')

global figK nu design_vec_iteration
figK = 0;
nu = 0.65; % reduced volume

%%% write results in the diary file
diary_name = ['max_efficiency_diary_nu_' num2str(nu*100,'%.3i') '.txt'];
if isfile(diary_name)
    delete(diary_name);
end
diary(diary_name);
fprintf('Date and Time: %s \n', datetime('now'));
fprintf('MATLAB Version: %s \n', version);
fprintf('Computer Type: %s \n', computer);
[~, chipInfo] = system('sysctl -n machdep.cpu.brand_string');
disp(['Chip Type: ', chipInfo]);

%%% START MAIN PROGRAM %%%
fprintf('-------START-------\n-------------------\n')
tStart = tic;

%%% Read Initial Shape from MinDragForce %%%
dir_name = '../min_drag_force/';
file1 = [dir_name 'final_designvec_mindragforce_nu_' num2str(100*nu, '%.3i') '.txt'];
file2 = [dir_name 'all_designvec_mindragforce_nu_' num2str(100*nu, '%.3i') '.txt'];
fID1 = fopen(file1, 'r'); 
fID2 = fopen(file2, 'r');
design_vec_temp = fscanf(fID1, '%f');
all_design_vec = fscanf(fID2, '%f');
dimvec = length(design_vec_temp);
fclose(fID1); fclose(fID2);
% goback = 1; % use the last result
% goback = 2; % use the 2nd last result
goback = 4;
design_vec = all_design_vec( end-goback*dimvec+1 : end-(goback-1)*dimvec );

design_vec_iteration = design_vec;
shape_initial = shape3Dmaxefficiency2(design_vec);
fprintf('Error between JD/JW and E, %.2e \n\n',shape_initial.JD/shape_initial.JW - shape_initial.JE)
plot(shape_initial.arclen,shape_initial.uslip)
saveas(gcf, ['./maxE_iteration_0_nu_' num2str(100*nu, '%.3i') '_uslip.png']);
close all

THETA_bas = matrix_theta_basis(shape_initial.t, shape_initial.NL, shape_initial.L);
fprintf('\n-- Print Design Parameters: --\n')
fprintf('p:%i,np:%i,NL:%i,NLuslip:%i,L:%g,Luslip:%g,Nu:%i\n',...
    shape_initial.p,shape_initial.np,shape_initial.NL,...
    shape_initial.NLuslip,shape_initial.L,shape_initial.Luslip,shape_initial.Nu);
shape_initial.plotblack;
saveas(gcf, ['./maxE_iteration_0_nu_' num2str(100*nu, '%.3i') '.png']);
close all

fprintf('\n-- Print Initial Results: --\n')
shape_initial.printresults;

%%% fixate poles
global fixed_RN fixed_RS fixed_ZN fixed_ZS
fixed_dim = 0; % if set to 0, not fixate the poles manually
[new_design_vec,fixed_RN,fixed_RS,fixed_ZN,fixed_ZS] = fixpoles(design_vec,fixed_dim);

%%% set optimization parameters
c.lam = 0; 
fprintf('--> c.lam: %f \n\n', c.lam);
c.sig = 15000;
fprintf('--> c.sig: %g \n\n', c.sig);
c.target = nu;
fprintf('--> c.target (reduced volume): %g \n\n', c.target);
c.tolerance = c.sig^(-0.1);
fprintf('--> c.tolerance: %g \n\n',c.tolerance);
c.multi_cst = 1;
fprintf('--> c.multi_cst: %g \n\n',c.multi_cst);
constraint_tolerance_nomult = 0.001;
fprintf('--> constraint_tolerance_nomult: %g \n\n',constraint_tolerance_nomult);
constraint_tolerance = constraint_tolerance_nomult*c.multi_cst;
fprintf('--> constraint_tolerance: %g \n\n',constraint_tolerance);
increaseSIGfactor = 10;
fprintf('--> increaseSIGfactor: %g \n\n',increaseSIGfactor);

options = optimoptions(@fminunc, ...
    'Display','iter', ...
    'Algorithm','quasi-newton', ...
    'HessianApproximation','lbfgs', ...
    'SpecifyObjectiveGradient',true, ...
    'UseParallel',true, ...
    'StepTolerance', 1e-3, ...
    'OutputFcn',@outfun_maxE ...
    );

L_A = @(vec) objective_maxE(vec,THETA_bas,c,'grad on',fixed_RN,fixed_RS,fixed_ZN,fixed_ZS);
[new_design_vec_update,fval,exitflag,output] = fminunc(L_A,new_design_vec,options);

new_design_vec = new_design_vec_update;
design_vec = fixpolesreverse(new_design_vec,fixed_RN,fixed_RS,fixed_ZN,fixed_ZS);

shape_current = shape3Dmaxefficiency2(design_vec);

fprintf('Error between JD/JW and E, %.2e \n\n',shape_current.JD/shape_current.JW - shape_current.JE)
plot(shape_current.arclen,shape_current.uslip)
saveas(gcf, './iteration_final_uslip.png');
close all

tEnd = toc(tStart);
fprintf('Total Elapsed Time: %i hours, %i minutes, %i seconds. \n\n', ...
    round(tEnd/3600), round(mod(tEnd,3600)/60), ceil(mod(tEnd,60)));

diary off