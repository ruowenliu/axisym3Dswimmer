% Read design vector and plot shape save in eps

close all
clear

%% nu=0.50
tic
file_name = '/Users/ruowen/Documents/GitHub_Ruowen/axisym3Dswimmer/min_drag_1/final_designvec_mindragforce_nu_5.00e-01.txt';
fID = fopen(file_name, 'r');
design_vec = fscanf(fID, '%f');
fclose(fID);
shape = shape3Dmaxefficiency2(design_vec);
fprintf('Verify U is one: %g\n', shape.U);
shape.plotgreen;
saveas(gcf,'mindrag_nu_050_opt_shape','epsc');
toc

%% nu=0.60
tic
file_name = '/Users/ruowen/Documents/GitHub_Ruowen/axisym3Dswimmer/min_drag_2/final_designvec_mindragforce_nu_060.txt';
fID = fopen(file_name, 'r');
design_vec = fscanf(fID, '%f');
fclose(fID);
shape = shape3Dmaxefficiency2(design_vec);
fprintf('Verify U is one: %g\n', shape.U);
shape.plotgreen;
saveas(gcf,'mindrag_nu_060_opt_shape','epsc');
toc

%% nu=0.70
tic
file_name = '/Users/ruowen/Documents/GitHub_Ruowen/axisym3Dswimmer/min_drag_2/final_designvec_mindragforce_nu_070.txt';
fID = fopen(file_name, 'r');
design_vec = fscanf(fID, '%f');
fclose(fID);
shape = shape3Dmaxefficiency2(design_vec);
fprintf('Verify U is one: %g\n', shape.U);
shape.plotgreen;
saveas(gcf,'mindrag_nu_070_opt_shape','epsc');
toc

%% nu=0.80
tic
file_name = '/Users/ruowen/Documents/GitHub_Ruowen/axisym3Dswimmer/min_drag_2/final_designvec_mindragforce_nu_080.txt';
fID = fopen(file_name, 'r');
design_vec = fscanf(fID, '%f');
fclose(fID);
shape = shape3Dmaxefficiency2(design_vec);
fprintf('Verify U is one: %g\n', shape.U);
shape.plotgreen;
saveas(gcf,'mindrag_nu_080_opt_shape','epsc');
toc

%% nu=0.90
tic
file_name = '/Users/ruowen/Documents/GitHub_Ruowen/axisym3Dswimmer/min_drag_2/final_designvec_mindragforce_nu_090.txt';
fID = fopen(file_name, 'r');
design_vec = fscanf(fID, '%f');
fclose(fID);
shape = shape3Dmaxefficiency2(design_vec);
fprintf('Verify U is one: %g\n', shape.U);
shape.plotgreen;
saveas(gcf,'mindrag_nu_090_opt_shape','epsc');
toc