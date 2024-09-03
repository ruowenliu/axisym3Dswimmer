% This file produces the shape optimization process from a prolate spheroid
% to the minimum drag force shape,
% constraned with different reduced volumes nu
% Written by Ruowen Liu, May 2023
% Last modified and tested in September 2024

%%% preset
clear; close all; format long; format compact;
warningid = 'MATLAB:nearlySingularMatrix';
warning('off',warningid)
addpath('../shape_classes')
addpath('../quadrature_and_kernal')

global nu

input_data = readtable('input_parameters.csv');

for row = size(input_data,1)

    nu = input_data{row, 1};

    %%% write results in the diary file
    diary_name = ['min_drag_force_diary_nu_' num2str(100*nu, '%.3i') '.txt'];
    if isfile(diary_name)
        delete(diary_name);
    end
    diary(diary_name);
    fprintf('Date and Time: %s \n', datetime('now'));
    fprintf('MATLAB Version: %s \n', version);
    fprintf('Computer Type: %s \n', computer);
    [~, chipInfo] = system('sysctl -n machdep.cpu.brand_string');
    disp(['Chip Type: ', chipInfo]);

    %%% the initial prolate spheroid
    rng('default');
    design_vec = get_initial_prolate(nu + 0.01*randn(1),'no noise');
    shape_initial = shape3Dmaxefficiency2(design_vec);
    shape_initial.printresults;
    fprintf('\nError between JD/JW and E, %.1e (It should be close to zero.)\n\n',shape_initial.JD/shape_initial.JW - shape_initial.JE)
    plot(shape_initial.arclen, shape_initial.uslip)
    saveas(gcf, ['./iteration_0_nu_' num2str(100*nu, '%.3i') '_uslip.png']);
    close all

    fprintf('-------START-------\n-------------------\n')

    tStart = tic;

    THETA_bas = matrix_theta_basis(shape_initial.t, shape_initial.NL, shape_initial.L);
    fprintf('--> Print Design Parameters: \n')
    fprintf('p:%i,np:%i,NL:%i,NLuslip:%i,L:%g,Luslip:%g,Nu:%i \n\n',shape_initial.p,shape_initial.np,shape_initial.NL,shape_initial.NLuslip,shape_initial.L,shape_initial.Luslip,shape_initial.Nu);

    c.lam = input_data{row, 2}; c.sig = input_data{row, 3}; 
    fprintf('--> c.lam: %.5f and c.sig: %.5f \n\n', c.lam, c.sig);

    c.target = nu; 
    fprintf('--> c.target (reduced volume): %.2f \n\n', c.target);

    c.tolerance = c.sig^(-0.1); 
    fprintf('--> c.tolerance: %g \n\n',c.tolerance);

    c.multi_cst = input_data{row, 4}; 
    fprintf('--> c.multi_cst: %g \n\n',c.multi_cst);

    constraint_tolerance = 1e-4 * c.multi_cst;
    fprintf('--> constraint_tolerance: %g \n\n',constraint_tolerance);

    increaseSIGfactor = input_data{row, 5}; 
    fprintf('--> increaseSIGfactor: %g \n\n',increaseSIGfactor);

    step_tolerance = input_data{row, 6}; 
    fprintf('--> StepTolerance: %g \n\n',step_tolerance);

    optimality_tolerance = input_data{row, 7}; 
    fprintf('--> OptimalityTolerance: %g \n\n',optimality_tolerance);
    
    options = optimoptions(@fminunc, ...
        'Display','iter', ...
        'Algorithm','quasi-newton', ...
        'HessianApproximation','lbfgs', ...
        'SpecifyObjectiveGradient',true, ...
        'UseParallel',true, ...
        'StepTolerance', step_tolerance, ...
        'OptimalityTolerance', optimality_tolerance, ...
        'OutputFcn',@outfunDragForceMin ...
        );

    stop_code = 0; % flag to stop optimization
    alm_num = -1; % counter
    outputiternum = 0; % counter

    while stop_code == 0

        alm_num = alm_num + 1; fprintf('\n**** ALM Loop = %i \n', alm_num);
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
    fprintf('Total Elapsed Time: %i hours, %i minutes, %i seconds. \n\n', ...
        round(tEnd/3600), round(mod(tEnd,3600)/60), ceil(mod(tEnd,60)));

    fprintf('Final Shape:\n')
    shape_current = shape3Dmaxefficiency2(design_vec);
    shape_current.printresults;
    fprintf('\nError between JD/JW and E, %.1e (It should be close to zero.)\n\n',shape_current.JD/shape_current.JW - shape_current.JE)
    plot(shape_current.arclen, shape_current.uslip)
    saveas(gcf, ['./iteration_final_nu_' num2str(100*nu, '%.3i') '_uslip.png']);
    close all

    fID = fopen(['./final_designvec_mindragforce_nu_' num2str(100*nu, '%.3i') '.txt'], 'w');
    fprintf(fID, '%.15f \n', design_vec);
    fclose(fID);

    diary off

end