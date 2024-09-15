% Read design vector and plot shape save in eps

%% nu=0.50
tic
file_name = '/Users/ruowen/Documents/GitHub_Ruowen/axisym3Dswimmer/maxE_nu_050etc/maxE_result_5.00e-01/maxE_nu_5.00e-01_designvec_iter_2.txt';
fID = fopen(file_name, 'r');
design_vec = fscanf(fID, '%f');
fclose(fID);
shape = shape3Dmaxefficiency2(design_vec);
fprintf('Verify U is one: %g\n', shape.U);
shape.plotorange;
saveas(gcf,'maxE_nu_050_opt_shape','epsc');
toc

%% nu=0.60
tic
file_name = '/Users/ruowen/Documents/GitHub_Ruowen/axisym3Dswimmer/maxE_nu_060/back2sig18000okay/maxE_nu_060_designvec_iteration_4.txt';
fID = fopen(file_name, 'r');
design_vec = fscanf(fID, '%f');
fclose(fID);
shape = shape3Dmaxefficiency2(design_vec);
fprintf('Verify U is one: %g\n', shape.U);
shape.plotorange;
saveas(gcf,'maxE_nu_060_opt_shape','epsc');
toc

%% nu=0.70
tic
file_name = '/Users/ruowen/Documents/GitHub_Ruowen/axisym3Dswimmer/maxE_from_peanut/design_vec_final.txt';
fID = fopen(file_name, 'r');
design_vec = fscanf(fID, '%f');
fclose(fID);
shape = shape3Dmaxefficiency2(design_vec);
fprintf('Verify U is one: %g\n', shape.U);
shape.plotorange;
saveas(gcf,'maxE_nu_070_opt_shape','epsc');
toc

%% nu=0.80
tic
file_name = '/Users/ruowen/Documents/GitHub_Ruowen/axisym3Dswimmer/maxE_nu_080/maxE_nu_8.00e-01_designvec_iter_2.txt';
fID = fopen(file_name, 'r');
design_vec = fscanf(fID, '%f');
fclose(fID);
shape = shape3Dmaxefficiency2(design_vec);
fprintf('Verify U is one: %g\n', shape.U);
shape.plotorange;
saveas(gcf,'maxE_nu_080_opt_shape','epsc');
toc

%% nu=0.90
tic
file_name = '/Users/ruowen/Documents/GitHub_Ruowen/axisym3Dswimmer/maxE_nu_090/maxE_nu_9.00e-01_designvec_iter_10.txt';
fID = fopen(file_name, 'r');
design_vec = fscanf(fID, '%f');
fclose(fID);
shape = shape3Dmaxefficiency2(design_vec);
fprintf('Verify U is one: %g\n', shape.U);
shape.plotorange;
saveas(gcf,'maxE_nu_090_opt_shape','epsc');
toc