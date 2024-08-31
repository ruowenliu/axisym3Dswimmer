% This file verifies the shape sensitivities,
% compared with central finite difference
% Written by Ruowen Liu, June 2023
% Last modified and tested in August 2024

%%% preset
clear; close all;
format long; format compact;
warningid = 'MATLAB:nearlySingularMatrix';
warning('off',warningid)
addpath('../shape_classes')
addpath('../quadrature_and_kernal')
%%% write results in the diary file
diary_name = 'sens_verify_result.txt';
if isfile(diary_name)
    delete(diary_name);
end
diary(diary_name);
fprintf('Date and Time: %s \n', datetime('now'));
fprintf('MATLAB Version: %s \n', version);
fprintf('Computer Type: %s \n', computer);
[~, chipInfo] = system('sysctl -n machdep.cpu.brand_string');
disp(['Chip Type: ', chipInfo]);
%%%
shapelist = ["long-lump","asymmetric-vase","three-bumps-2"];
nu = 0.7; % set nu only for "prolate"
type_0 = "symmetric-peanut";
hs = 1e-3; % stepsize for central finite difference scheme
fprintf('\nFinite Difference step size: %g \n\n', hs);
for t1 = 1:length(shapelist)
    tic    
    type_1 = shapelist(t1);
    if type_0~=type_1
        %%% original shape
        design_vec_0 = get_initial_parameterized_shape_gallery(type_0,nu);
        shape_0 = shape3Dmaxefficiency2(design_vec_0); % must find optimal uslip here
        shape_0.plotblack;
        %%% purturbed to new shape
        design_vec_1 = get_initial_parameterized_shape_gallery(type_1,nu);
        shape_1 = shape3Dmaxefficiency2(design_vec_1); % must find optimal uslip here
        shape_1.plotgreen;
        %%%
        delta_design_vec = design_vec_1 - design_vec_0;
        %%% upstream and downstream for central finite difference scheme
        design_vec_upstream = design_vec_0 - hs*delta_design_vec;
        design_vec_downstream = design_vec_0 + hs*delta_design_vec;
        shape_upstream = shape3Dmaxefficiency2(design_vec_upstream);
        shape_downstream = shape3Dmaxefficiency2(design_vec_downstream);
        %%%
        tht = shape_1.x - shape_0.x;
        shapeDer = shape_derivatives(shape_0, tht);
        checklist = {'Jdrag_rByV','JE'};
        %%%
        fprintf('-----------------------------------------\n');
        fprintf("-- from %s to %s --\n", type_0, type_1);
        fprintf('-----------------------------------------\n');
        for k = 1:length(checklist)
            if strcmp(checklist{k}, 'Jdrag_rByV')
                rename = 'Jdrag';
            elseif strcmp(checklist{k}, 'JE')
                rename = 'E';
            else
                rename = checklist{k};
            end
            quantity_name = checklist{k};
            quantity_name_prime = [quantity_name 'prime'];
            prime_FD = (shape_downstream.(quantity_name)-shape_upstream.(quantity_name))/2/hs;
            fprintf('%s of initial shape: %.10f \n\n', rename, shape_0.(quantity_name));
            fprintf('check derivative %s''\n\n', rename);
            fprintf('The derivative %s'' by central finite difference:   %.10f \n', rename, prime_FD);
            fprintf('The derivative %s'' by derived sensitivity formula: %.10f \n\n', rename, shapeDer.(quantity_name_prime));
            fprintf('Absolute error of %s'': %.2e \n', ...
                rename, abs(prime_FD-shapeDer.(quantity_name_prime)));
            fprintf('Relative error of %s'': %.2e \n\n', ...
                rename, abs((prime_FD-shapeDer.(quantity_name_prime))/prime_FD));
        end
        elapsedTime = toc;  % End the timer and capture the elapsed time
        fprintf('Elapsed time for one comparison: %.1f seconds \n\n', elapsedTime);
    end
end

diary off