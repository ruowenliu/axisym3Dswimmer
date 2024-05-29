% (c) 05/2023 Ruowen Liu
% Just compute properties of prolate spherioid with given reduced volume nu

clear; close all; format long; format compact;
warningid = 'MATLAB:nearlySingularMatrix';
warning('off',warningid)

for nu = 0.6:0.05:1.0 % set reduced volume

    tic
    disp(nu)
    ProlateSpheroidAxisymSwimmer(nu)

    toc
end


% %%
% nu = 0.6:0.05:1.0;
% for imk = 1:9
%     nu_im = nu(imk);
%     h = figure(imk);
%     saveas(h,['ProlateSpheroidData' num2str(nu_im*100,'%2i')],'epsc');
% end
% 
