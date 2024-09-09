%%% Update c
% Cnu = c.multi_cst*(shape_current.rvol - c.target);
%
% c.lam = c.lam - c.sig*(Cnu);
% c.tolerance = c.sig^(-0.9)*c.tolerance;
%
% or
c.sig = increaseSIGfactor*c.sig;
% c.tolerance = c.sig^(-0.1);

% 
% disp(c);
% keyboard

%%% START MAIN PROGRAM %%%
fprintf('-------START-------\n-------------------\n')
tStart = tic;

%%% Read Initial Shape from MinDragForce %%%
file1 = './back3sig100/maxE_nu_065_designvec_iteration_3.txt';
% file2 = [dir_name 'all_designvec_mindragforce_nu_' num2str(100*nu, '%.3i') '.txt'];
fID1 = fopen(file1, 'r'); 
% fID2 = fopen(file2, 'r');
design_vec = fscanf(fID1, '%f');
% all_design_vec = fscanf(fID2, '%f');
% dimvec = length(design_vec_temp);
fclose(fID1); 
% fclose(fID2);
% % goback = 1; % use the last result
% % goback = 2; % use the 2nd last result
% goback = 3; % use the 3rd last result
% design_vec = all_design_vec( end-goback*dimvec+1 : end-(goback-1)*dimvec );

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
% c.lam = ;
% c.sig = ;
fprintf('--> c.lam: %.4f and c.sig: %g \n\n', c.lam, c.sig);
c.target = nu;
fprintf('--> c.target (reduced volume): %g \n\n', c.target);
% c.tolerance = ;
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
    'StepTolerance', 1e-4, ...
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