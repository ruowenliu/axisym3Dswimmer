% This file produces the shape optimization process from a prolate spheroid
% to the minimum drag force shape,
% constraned with reduced volume nu=0.60
% Written by Ruowen Liu, May 2023
% Last modified and tested in September 2024

%%% preset
clear; close all; format long; format compact;
warningid = 'MATLAB:nearlySingularMatrix';
warning('off',warningid)
addpath('../shape_classes')
addpath('../quadrature_and_kernal')

%%% the initial prolate spheroid
nu = 0.60;
rng('default');
design_vec = get_initial_prolate(nu + 0.01*randn(1),'no noise');
shape_initial = shape3Dadjoint(design_vec);
shape_initial.printresults;

%%% write results in the diary file
diary_name = 'min_drag_main_nu_60e-2.txt';
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

tStart = tic;

THETA_bas = matrix_theta_basis(shape_initial.t, shape_initial.NL, shape_initial.L);
fprintf('--> Print Design Parameters: \n')
fprintf('p:%i,np:%i,NL:%i,NLuslip:%i,L:%g,Luslip:%g,Nu:%i \n\n',shape_initial.p,shape_initial.np,shape_initial.NL,shape_initial.NLuslip,shape_initial.L,shape_initial.Luslip,shape_initial.Nu);

c.lam = 0; c.sig = 8; fprintf('--> c.lam: %.5f and c.sig: %.5f \n\n', c.lam, c.sig);

c.target = nu; fprintf('--> c.target (reduced volume): %.2f \n\n', c.target);

c.tolerance = c.sig^(-0.1); fprintf('--> c.tolerance: %g \n\n',c.tolerance);

c.multi_cst = 2; fprintf('--> c.multi_cst: %g \n\n',c.multi_cst);

constraint_tolerance = 1e-3*2^(-2); fprintf('--> constraint_tolerance: %g \n\n',constraint_tolerance);

increaseSIGfactor = 1.5; fprintf('--> increaseSIGfactor: %g \n\n',increaseSIGfactor);

options = optimoptions(@fminunc, ...
    'Display','iter', ...
    'Algorithm','quasi-newton', ...
    'HessianApproximation','lbfgs', ...
    'SpecifyObjectiveGradient',true, ...
    'UseParallel',true, ...
    'StepTolerance', 1e-4, ...
    'OptimalityTolerance', 1e-4, ...
    'OutputFcn',@outfunDragForceMin ...
    ); 

stop_code = 0; % flag to stop optimization
alm_num = -1; % counter
outputiternum = 0; % counter

while stop_code == 0

    alm_num = alm_num + 1; fprintf('**** ALM Loop = %i \n', alm_num);
    fprintf(['current parameters: lambda = ', num2str(c.lam,4), ' and sigma = ' num2str(c.sig), '\n']);
    fprintf('current c.tolerance: %g \n\n', c.tolerance);

    L_A = @(vec) objectiveDragForceMin(vec, THETA_bas, c, 'grad on');
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
        c.lam = 0;
        c.tolerance = c.sig^(-0.1);
    end

    if isempty(output.stepsize)
        outputiternum = outputiternum + 1;
    end

    if outputiternum >= 10
        stop_code = 1; 
        fprintf('-----> Optimization Halted \n\n')
        
    end
end

tEnd = toc(tStart);
fprintf('Total Elapsed Time: %i hours, %i minutes, %i seconds. \n\n', ...
    round(tEnd/3600), round(mod(tEnd,3600)/60), ceil(mod(tEnd,60)));

shape_current = shape3Dmaxefficiency(design_vec);
shape_current.printresults;

fprintf('Final Shape:\n')
shape_current = shape3Dmaxefficiency2(design_vec);
shape_current.printresults;

diary off