function ProlateSpheroidAxisymSwimmer(nu)
% (c) 05/2023 Ruowen Liu

fprintf('---------------------\n---------------------\n')

diaryname = ['ProlateSpheroidDiary' num2str(nu*100,'%2i') '.txt'];
diary(diaryname);
diary on

fprintf('  Calculate Prolate Spheroid Shape by Bslpines.\n')
design_vec = get_initial_prolate(nu,'no noise'); % or 'add noise'
shape_initial = shape3Dmaxefficiency(design_vec); % must find optimal uslip here

fprintf('    Print Shape Design Parameters: --\n')
fprintf('p:%i,np:%i,NL:%i,NLuslip:%i,L:%g,Luslip:%g,Nu:%i\n',shape_initial.p,shape_initial.np,shape_initial.NL,shape_initial.NLuslip,shape_initial.L,shape_initial.Luslip,shape_initial.Nu);

fprintf('    Results:\n')
shape_initial.printresults;
shape_initial.plotblue;
saveas(gcf,['ProlateSpheroidData' num2str(nu*100,'%2i')],'pdf')

savedatafilename = ['ProlateSpheroidData' num2str(nu*100,'%2i') '.mat'];
save(savedatafilename)

savefigfilename = ['ProlateSpheroidFigure' num2str(nu*100,'%2i')];
saveas(gcf, savefigfilename, 'fig')

diary off
