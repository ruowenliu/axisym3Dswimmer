% This file produces the shape optimization process from a peanut-like
% shape to the maximum swimming efficiency shape,
% constraned with reduced volume nu=0.7
% Written by Ruowen Liu, July 2023
% Last modified and tested in August 2024

%%% preset
clear; close all; format long; format compact;
warningid = 'MATLAB:nearlySingularMatrix';
warning('off',warningid)
addpath('../shape_classes')
addpath('../quadrature_and_kernal')

global figK nu design_vec_iteration
figK = 10; %%%%%%%%
nu = 0.65; % reduced volume

%%% write results in the diary file
diary_name = ['max_efficiency_nu_' num2str(nu*100,'%.3i') '.txt'];
% if isfile(diary_name)
%     delete(diary_name);
% end
diary(diary_name);
fprintf('Date and Time: %s \n', datetime('now'));
fprintf('MATLAB Version: %s \n', version);
fprintf('Computer Type: %s \n', computer);
[~, chipInfo] = system('sysctl -n machdep.cpu.brand_string');
disp(['Chip Type: ', chipInfo]);

%%% START MAIN PROGRAM %%%
fprintf('-------START-------\n-------------------\n')
tStart = tic;

fprintf('Calculate Initial Shape.\n')
% design_vec = get_initial_parameterized_shape_gallery("symmetric-peanut",nu);

% fileID = fopen(['initial' num2str(nu*100,'%.3i') '.txt'], 'r');
fileID = fopen('designvec_iteration_5.txt', 'r');
design_vec = fscanf(fileID, '%f');
fclose(fileID);
fprintf('\n\n Initial Design Vector: \n')
disp(design_vec);

design_vec_iteration = design_vec;
shape_initial = shape3Dmaxefficiency2(design_vec);
disp(shape_initial.JD/shape_initial.JW - shape_initial.JE)
plot(shape_initial.arclen,shape_initial.uslip)
keyboard


THETA_bas = matrix_theta_basis(shape_initial.t, shape_initial.NL, shape_initial.L);
fprintf('\n-- Print Design Parameters: --\n')
fprintf('p:%i,np:%i,NL:%i,NLuslip:%i,L:%g,Luslip:%g,Nu:%i\n',...
    shape_initial.p,shape_initial.np,shape_initial.NL,...
    shape_initial.NLuslip,shape_initial.L,shape_initial.Luslip,shape_initial.Nu);
shape_initial.plotblack;
saveas(gcf, './iteration_0.png');
close all
fprintf('\n-- Print Initial Results: --\n')
shape_initial.printresults;

%%% fixate poles
global fixed_RN fixed_RS fixed_ZN fixed_ZS
fixed_dim = 4; [new_design_vec,fixed_RN,fixed_RS,fixed_ZN,fixed_ZS] = fixpoles(design_vec,fixed_dim);

%%% set optimization parameters
c.lam = 0; 
c.sig = 10000; %c.sig = 200; 
fprintf('--> c.lam: %.4f and c.sig: %g \n\n', c.lam, c.sig);
c.target = nu; 
fprintf('--> c.target (reduced volume): %g \n\n', c.target);
c.tolerance = 0.1* c.sig^(-0.1); 
fprintf('--> c.tolerance: %g \n\n',c.tolerance);
c.multi_cst = 1; 
fprintf('--> c.multi_cst: %g \n\n',c.multi_cst);
constraint_tolerance_nomult = 0.0005; 
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

stop_code = 0; % flag to stop optimization
alm_num = -1; % counter
outputiternum = 0; % counter

while stop_code == 0
    alm_num = alm_num + 1; fprintf('\n ------ ALM Loop = %i \n\n', alm_num);

    % switch alm_num
    %     case {0,1}
    %         options.MaxIterations = 4;
    %         fprintf('**** set options.MaxIterations: %g\n',options.MaxIterations);
    %     case {2,3,4}
    %         options.MaxIterations = 5;
    %         fprintf('**** set options.MaxIterations: %g\n',options.MaxIterations);
    %         fixed_dim = 0;[new_design_vec,fixed_RN,fixed_RS,fixed_ZN,fixed_ZS] = fixpoles(design_vec,fixed_dim);
    %     otherwise
            options.MaxIterations = 100;
            fprintf('**** set options.MaxIterations: %g\n',options.MaxIterations);
            % go back to non-fix top situation
            fixed_dim = 0;[new_design_vec,fixed_RN,fixed_RS,fixed_ZN,fixed_ZS] = fixpoles(design_vec,fixed_dim);
    % end

    fprintf(['* In this ALM, lambda = ', num2str(c.lam,4), ' and sigma = ' num2str(c.sig), '\n']);
    fprintf('** c.tolerance: %g \n',c.tolerance);
    fprintf('*** fixate poles basis number: %i \n', fixed_dim);


    L_A = @(vec) objective_maxE(vec,THETA_bas,c,'grad on',fixed_RN,fixed_RS,fixed_ZN,fixed_ZS);
    [new_design_vec_update,fval,exitflag,output] = fminunc(L_A,new_design_vec,options);

    fprintf(' ------ ALM Loop = %i finished. \n\n', alm_num);

    new_design_vec = new_design_vec_update;
    design_vec = fixpolesreverse(new_design_vec,fixed_RN,fixed_RS,fixed_ZN,fixed_ZS);

    shape_current = shape3Dbasic(design_vec);
    Cnu = c.multi_cst*(shape_current.rvol - c.target);

    % save(['./' resultdir '/data_ALM_done_' num2str(alm_num) '.mat'])

    if abs(Cnu) < c.tolerance
        if abs(Cnu) < constraint_tolerance
            fprintf('-----> Optimization is complete. \n\n')
            stop_code = 1;
        else
            c.lam = c.lam - c.sig*(Cnu);
            c.tolerance = 0.1* c.sig^(-0.9)*c.tolerance;
        end
    else
        c.sig = increaseSIGfactor*c.sig;
        c.tolerance = 0.1* c.sig^(-0.1);
    end

    if isempty(output.stepsize)
        outputiternum = outputiternum + 1;
    end

    if outputiternum >= 2
        stop_code = 1; 
        fprintf('----->  Optimization done. Please investigate the results. \n\n')
    end

end

tEnd = toc(tStart);
fprintf('Total Elapsed Time: %i hours, %i minutes, %i seconds. \n\n', ...
    round(tEnd/3600), round(mod(tEnd,3600)/60), ceil(mod(tEnd,60)));

fprintf('\n\n Final Design Vector: \n')
disp(design_vec)

diary off