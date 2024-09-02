% This file produces the shape optimization process 
% from the shape with minimum drag force of a constrained reduced volume
% to the maximum swimming efficiency
% Written by Ruowen Liu, May 2023
% Last modified and tested in September 2024

%%% preset
clear; close all; format long; format compact;
warningid = 'MATLAB:nearlySingularMatrix';
warning('off',warningid)
addpath('../shape_classes')
addpath('../quadrature_and_kernal')

%%% read the shape with minimum drag force
nu = 0.70;
dir_name = '../min_drag_force_various_nu/';
file_name = [dir_name 'min_drag_main_nu_' num2str(100*nu, '%.3i') '.mat'];
load(file_name, 'design_vec');
design_vec_0 = design_vec;
shape_initial = shape3Dmaxefficiency2(design_vec_0);
shape_initial.plotgreen;

%%% write results in the diary file
diary_name = ['max_efficiency_main_nu_' num2str(100*nu, '%.3i') '.txt'];
if isfile(diary_name)
    delete(diary_name);
end
diary(diary_name);
fprintf('Date and Time: %s \n', datetime('now'));
fprintf('MATLAB Version: %s \n', version);
fprintf('Computer Type: %s \n', computer);
[~, chipInfo] = system('sysctl -n machdep.cpu.brand_string');
disp(['Chip Type: ', chipInfo]);

fprintf('-------START-------\n-------------------\n')

%global design_vec_iterations result_iterations nu
%global design_vec_fixed_Rnorth design_vec_fixed_Rsouth design_vec_fixed_Znorth design_vec_fixed_Zsouth new_design_vec_half_dim

tStart = tic;

THETA_bas = matrix_theta_basis(shape_initial.t, shape_initial.NL, shape_initial.L);
fprintf('--> Print Design Parameters: \n')
fprintf('p:%i,np:%i,NL:%i,NLuslip:%i,L:%g,Luslip:%g,Nu:%i \n\n',shape_initial.p,shape_initial.np,shape_initial.NL,shape_initial.NLuslip,shape_initial.L,shape_initial.Luslip,shape_initial.Nu);

fprintf('--> Print Initial Results: \n\n')
shape_initial.printresults;

c.lam = 0; c.sig = 10000; 
fprintf('--> c.lam: %.5f and c.sig: %.5f \n\n', c.lam, c.sig);

c.target = nu; 
fprintf('--> c.target (reduced volume): %g \n\n', c.target);

c.tolerance = c.sig^(-0.1); 
fprintf('--> c.tolerance: %g \n\n',c.tolerance);

c.multi_cst = 1; 
fprintf('--> c.multi_cst: %g \n\n',c.multi_cst);

constraint_tolerance = 1e-4 * c.multi_cst; 
fprintf('--> constraint_tolerance: %g \n\n',constraint_tolerance);

increaseSIGfactor = 10; 
fprintf('--> increaseSIGfactor: %g \n\n',increaseSIGfactor);

options = optimoptions(@fminunc, ...
    'Display','iter', ...
    'Algorithm','quasi-newton', ...
    'HessianApproximation','lbfgs', ...
    'SpecifyObjectiveGradient',true, ...
    'UseParallel',true);

stop_code = 0; % flag to stop optimization
alm_num = -1; % counter
outputiternum = 0; % counter

while stop_code == 0
    
    alm_num = alm_num + 1; fprintf('\n**** ALM Loop = %i \n', alm_num);
        fprintf(['current parameters: lambda = ', num2str(c.lam,4), ' and sigma = ' num2str(c.sig), '\n']);
        fprintf('current c.tolerance: %g \n\n', c.tolerance);

    L_A = @(vec) objectiveMaxEffi(vec, THETA_bas, c, 'grad on');
    [design_vec_update,~,~,output] = fminunc(L_A,design_vec,options);

    design_vec = design_vec_update;

    shape_current = shape3Dbasic(design_vec);

    Cnu = c.multi_cst*(shape_current.rvol - c.target);
        if abs(Cnu) < c.tolerance
            if abs(Cnu) < constraint_tolerance
                fprintf('-----> Optimization is complete. \n\n')
                stop_code = 1;
            else
                c.lam = c.lam - c.sig*(Cnu);
                c.tolerance = c.sig^(-0.9)*c.tolerance;
            end
        else
            c.sig = increaseSIGfactor*c.sig;
            c.tolerance = c.sig^(-0.1);
        end

        if isempty(output.stepsize)
            outputiternum = outputiternum + 1;
        end

        if outputiternum >= 3
            stop_code = 1;
            fprintf('-----> Optimization Halted \n\n')
        end

end

tEnd = toc(tStart);
fprintf('Total Elapsed Time: %i hours, %i minutes, %i seconds. \n', ...
    round(tEnd/3600), round(mod(tEnd,3600)/60), ceil(mod(tEnd,60)));

fprintf('Final Shape:\n')
shape_current = shape3Dmaxefficiency2(design_vec);
shape_current.printresults;
shape_current.plotorange;

fprintf('\n\n Final Design Vector: \n')
disp(design_vec)

diary off


% 
% 
% tStart = tic;
% 
% c.target = nu;
% c.lam = 0;
% fprintf('* Reset c.sig %g \n',c.sig);
% c.multi_cst = 1;
% 
% filelocation=['/Users/ruowen/Documents/MATLAB/shapeoptswimmer2023/swimmer-con/DragForceMin/' dirname '/finaldata.mat'];
% load(filelocation,'design_vec')
% design_vec_0 = design_vec;
% 
% fprintf('--- Initial ---\n');
% shape_initial = shape3Dmaxefficiency(design_vec_0);
% shape_initial.printresults;
% shape_initial.plotgreen;
% saveas(gcf, ['./' resultdir '/shape_'  num2str(nu*100,'%2i') '_initial_dragmin'],'pdf')
% saveas(gcf, ['./' resultdir '/shape_'  num2str(nu*100,'%2i') '_initial_dragmin'],'epsc')
% THETA_bas = matrix_theta_basis(shape_initial.t, shape_initial.NL, shape_initial.L);
% 
% options = optimoptions(@fminunc, ...
%     'Display','iter', ...
%     'Algorithm','quasi-newton', ...
%     'HessianApproximation','lbfgs', ...
%     'SpecifyObjectiveGradient',true, ...
%     'UseParallel',true ...    
%     ); 
% %    'OutputFcn',@outfunMaxEffi, ...
% % 'StepTolerance', 1e-3 ...
% 
% L_A = @(vec) objectiveMaxEffi(vec, THETA_bas, c, 'grad on');
% [design_vec_final,~,~,output] = fminunc(L_A,design_vec_0,options);
% 
% shape_final = shape3Dmaxefficiency(design_vec_final);
% shape_final.printresults;
% shape_final.plotorange;
% 
% save(['./' resultdir '/data_' num2str(nu*100,'%2i') '_final.mat']);
% 
% saveas(gcf, ['./' resultdir '/shape_'  num2str(nu*100,'%2i') '_final_maxeffi'],'pdf')
% saveas(gcf, ['./' resultdir '/shape_'  num2str(nu*100,'%2i') '_final_maxeffi'],'epsc')
% 
% tEnd = toc(tStart);
% fprintf('Total Elapsed Time: %i hours, %i minutes, %i seconds. \n', ...
%     round(tEnd/3600), round(mod(tEnd,3600)/60), ceil(mod(tEnd,60)));
% 
% close all
% diary off
% 
% end