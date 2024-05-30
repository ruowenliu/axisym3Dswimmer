function design_vec = get_initial_prolate(nu,noise)
% (c) 05/30/2024 Ruowen Liu

if nu==1
    axis.a = 1; axis.b = 1; % sphere
else
    find_ab = @(ab) [volume_prolate(ab(1),ab(2))-1,6*sqrt(pi)*volume_prolate(ab(1),ab(2))/(surface_area_prolate(ab(1),ab(2))^1.5)-nu];
    ab_val = fsolve(find_ab, [0.5,1]);
    axis.a = ab_val(1); axis.b = ab_val(2); % semi-major axis.a; semi-minor axis.b
end
shape_prolate = shape3Dparam; % Use np and NL in "shape3Dparam"
[xi_R_2L, xi_Z_2L] = prolate_Bsp(shape_prolate.np, axis, shape_prolate.NL);
xi_R = xi_R_2L(1:shape_prolate.NL+5); xi_Z = xi_Z_2L(1:shape_prolate.NL+5);
design_vec = [xi_R(2:end-1); xi_Z(2:end-1)]; % design vector for shape

% Add white noise to the design_vector
noise_level = 50; % the larger number the lower noise
if strcmp(noise,'add noise')
    design_vec = awgn(design_vec,noise_level); 
end

end
