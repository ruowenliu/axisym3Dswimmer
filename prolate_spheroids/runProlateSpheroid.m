% (c) 05/30/2024 Ruowen Liu
% Compute properties of prolate spherioid with given reduced volume nu

clear; close all; format long; format compact;
warningid = 'MATLAB:nearlySingularMatrix';
warning('off',warningid)

nu = 0.6; % or use a for loop for nu = 0.6:0.05:1.0
tic
fprintf('Prolate Spheroid, Reduced Volume = %g \n', nu)
ProlateSpheroidAxisymSwimmer(nu)
fprintf('Time elapsed: %g seconds.\n', toc)