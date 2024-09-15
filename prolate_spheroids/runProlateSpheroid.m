% (c) 05/30/2024 Ruowen Liu
% Compute properties of prolate spherioid with given reduced volume nu

clear; close all; format long; format compact;
warningid = 'MATLAB:nearlySingularMatrix';
warning('off',warningid)
addpath('../shape_classes')
addpath('../quadrature_and_kernal')

for nu = 0.50:0.05:1.00
    tic
    fprintf('Prolate Spheroid, Reduced Volume = %.2f \n', nu)

    fprintf('--- Calculate Prolate Spheroid Shape by Bslpines ---\n')
    design_vec = get_initial_prolate(nu,'no noise'); % or 'add noise' to shape
    shape_initial = shape3Dmaxefficiency2(design_vec); % find optimal uslip

    fprintf('--- Print Shape Design Parameters ---\n')
    fprintf('p:%i, np:%i, NL:%i, NLuslip:%i, L:%g, Luslip:%g, Nu:%i\n',...
        shape_initial.p, shape_initial.np, shape_initial.NL,...
        shape_initial.NLuslip, shape_initial.L, shape_initial.Luslip,...
        shape_initial.Nu);

    fprintf('--- Results ---\n')
    shape_initial.printresults;
    shape_initial.plotblue;

    savefigfilename = ['ProlateSpheroidFigure' num2str(round(nu*100),'%.3i')];
    saveas(gcf, savefigfilename, 'png')
    saveas(gcf, savefigfilename, 'epsc')

    disp(savefigfilename);
    shape_initial.printresults;

    fprintf('Time elapsed: %g seconds. \n\n', toc)

end