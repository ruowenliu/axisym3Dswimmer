% Copyright @ Ruowen Liu, July 2023
% Modified in Sep, 2024

close all
clear
warningid = 'MATLAB:nearlySingularMatrix'; warning('off',warningid);
addpath('../shape_classes')
addpath('../quadrature_and_kernal')

nu = 0.50;

for kcase = 1:3
    tic
    %%% load shape properties
    switch kcase

        case 3
            titleline = sprintf('Max Efficiency');

            mat_name = ['vfield_' num2str(kcase) '.mat'];
            if exist(mat_name, 'file') == 2
                disp('\n Load Results \n');
                load(mat_name);
            else
                disp('\n Recalculate Results \n');
                % Read design vector
                dir_result = ['../maxE_nu_050etc/maxE_result_' num2str(nu, '%.2e') '/'];
                file_name = [dir_result 'maxE_nu_5.00e-01_designvec_iter_2.txt'];
                fID = fopen(file_name, 'r');
                design_vec = fscanf(fID, '%f');
                fclose(fID);
                shape = shape3Dmaxefficiency2(design_vec);
                fprintf('Verify U is one: %g\n', shape.U);

            end

        case 2
            titleline = sprintf('Min Drag Force');
            mat_name = ['vfield_' num2str(kcase) '.mat'];
            if exist(mat_name, 'file') == 2
                disp('\n Load Results \n');
                load(mat_name);
            else
                disp('\n Recalculate Results \n');
                % Read design vector
                dir_result = ['../min_drag_1' '/'];
                file_name = [dir_result 'final_designvec_mindragforce_nu_5.00e-01.txt'];
                fID = fopen(file_name, 'r');
                design_vec = fscanf(fID, '%f');
                fclose(fID);
                shape = shape3Dmaxefficiency2(design_vec);
                fprintf('Verify U is one: %g\n', shape.U);
            end

        case 1
            titleline = sprintf('Prolate Spheroid');
            mat_name = ['vfield_' num2str(kcase) '.mat'];
            if exist(mat_name, 'file') == 2
                disp('\n Load Results \n');
                load(mat_name);
            else
                disp('\n Recalculate Results \n');
                % Read design vector
                design_vec = get_initial_prolate(nu,'no noise');
                shape = shape3Dmaxefficiency2(design_vec);
                fprintf('Verify U is one: %g\n', shape.U);
            end

    end

    fprintf('The north pole location: \n');
    disp(imag(shape.x(1)));

    if exist(mat_name, 'file') ~= 2

        %%% prepare plot 3D swimmer view
        xlimval = [-2,2]; ylimval = [-3,3]; zlimval = [-3.5,3.5];
        x = real(shape.x);
        y = imag(shape.x);
        xy = [x, y];
        r = xy(1:3:end,1); z = xy(1:3:end,2);
        n = 2000;
        t = linspace(0,2*pi,n);
        X = r*cos(t);
        Y = repmat(z,1,n);
        Z = r*sin(t);

        %%% Calculate velocity flow
        nr = 360; nz = nr/2*3;
        gr = linspace(xlimval(1), xlimval(2), nr+1);
        gz = linspace(ylimval(1), ylimval(2), nz+1); % notice using ylim (x->r, and y->z)
        [rr, zz] = meshgrid(gr,gz);
        xz = (rr+1i*zz);
        xpart = [real(shape.x); -real(flip(shape.x)); real(shape.x(1))];
        zpart = [imag(shape.x); imag(flip(shape.x)); imag(shape.x(1))];
        [IN, ~] = inpolygon(rr,zz,xpart,zpart); % if inside the swimmer, labeled 1

        % Manually exclude near-boundary domain %
        IN_new = IN;
        for rowk = 10:nz-10 % row corresponds to y (i.e. z)
            for colk = 10:nr-10 % column corresponds to x (i.e. r)
                boundval_lr = max(int8([IN(rowk,colk-1),IN(rowk,colk+1)]));
                boundval_ud = max(int8([IN(rowk-1,colk),IN(rowk+1,colk)]));
                %boundval_1 = max(int8([IN(rowk-1,colk),IN(rowk+1,colk),IN(rowk,colk-1),IN(rowk,colk+1),IN(rowk-1,colk-1),IN(rowk+1,colk+1),IN(rowk+1,colk-1),IN(rowk-1,colk+1)]));
                %boundval_2 = max(int8([IN(rowk-2,colk),IN(rowk-2,colk-1),IN(rowk-2,colk-2),IN(rowk-2,colk+1),IN(rowk-2,colk+2),...
                %             IN(rowk+2,colk),IN(rowk+2,colk-1),IN(rowk+2,colk-2),IN(rowk+2,colk+1),IN(rowk+2,colk+2),...
                %             IN(rowk,colk-2),IN(rowk-1,colk-2),IN(rowk+1,colk-2),IN(rowk,colk+2),IN(rowk-1,colk+2),IN(rowk+1,colk+2)]));
                boundval = boundval_lr+boundval_ud;
                if boundval > 0
                    if ~IN(rowk,colk)
                        IN_new(rowk,colk)=true;
                    end
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ii = ~IN_new; % if inside the swimmer, labeled 0; if outside the swimmer, labeled 1
        tstruct.x = xz(ii(:));  % stack ii to one vector (column 1; column 2; column 3; ...)
        U = shape.U; fprintf("Value of U: if the following number is 1, correct.\n"); disp(U);% must confirm that U=1
        rhs = [shape.uslip; shape.uslip].*[real(shape.tang);imag(shape.tang)]+[zeros(size(shape.t));U*ones(size(shape.t))];   % rhs
        A = AlpertSphereSLPMat(shape);
        mu =  A\rhs;                   % density
        u = nan*(1+1i)*xz;
        temp = SphereQuadtp(shape,tstruct)*mu;
        u(ii) = temp(1:end/2) + 1i*temp(end/2+1:end);  % evaluate velocity field

        % Plot velocity flow in lab frame
        type = 'lab_frame';
        u_magnitude = sqrt(real(u).^2 + imag(u).^2);

        max_mag = max(max(u_magnitude)); min_mag = min(min(u_magnitude));
        fprintf('\n %s min and max magnitude = [ %g, %g ] \n', type, min_mag, max_mag);

        save(['vfield_' num2str(kcase) '.mat']);
        timecompute = toc;
        fprintf('->  Generating data time cost about %.0f minutes. \n', round(timecompute/60)); % about 5 min each

    else

        %%% plotting
        tic
        % plot 3D swimmer view
        figure('Position', [0 0 600 720])
        hold on
        s = surf(X,Y,Z);
        s.EdgeAlpha = 0.0;
        s.FaceColor = 0.9*[1 1 1];
        s.FaceAlpha = 0.8;
        s.FaceLighting = 'gouraud';
        s.BackFaceLighting = 'unlit';
        axis equal
        grid off
        view(0,90)
        light('Position',[-0.4 0.2 0.3]); % Add lights
        xlim(xlimval); ylim(ylimval); zlim(zlimval);

        % plot stream
        hold on
        pc = pcolor(gr,gz,u_magnitude);
        set(pc, 'FaceColor', 'interp', 'LineStyle', 'none', 'FaceLighting', 'none')
        bm = brewermap([], '-Spectral');
        colormap(bm);
        clb = colorbar;cluplim=0.5;clim([0,cluplim]);clb.Ticks = 0:0.05:cluplim;clb.TickLabels{end}="$\ge "+num2str(cluplim)+"$";set(clb,'TickLabelInterpreter','latex');
        %if kcase < 3, clb.Visible = 'off'; end
        streamlinecolor = 'w'; streamlinewidth = 0.8;
        streamlinedensity = 1.8;
        % The default value is 1. Higher values produce more streamlines on each plane. For example, 2 produces approximately twice as many streamlines as the default, while 0.5 produces approximately half as many.
        h0 = streamslice(gr,gz,real(u),imag(u),streamlinedensity);
        set(h0,'color',streamlinecolor,'linewidth',streamlinewidth);
        h1 = streamslice(-gr,gz,-real(u),imag(u),streamlinedensity);
        set(h1,'color',streamlinecolor,'linewidth',streamlinewidth);

        % set tick and title
        titlefontsize = 40; tickfontsize = 35; xlabfontsize = 35;
        set(gca,'fontsize',tickfontsize,'TickLabelInterpreter','latex');
        set(gca,'XTick',xlimval(1):1:xlimval(end));
        set(gca,'YTick',ylimval(1):1:ylimval(end));
        box on

        lineA = ['$\nu$ ' num2str(shape.rvol,'%.3f')];lineB = ['$E$ ' num2str(shape.JE,'%.3f')];lineC = ['$J_{\! drag}$ ' num2str(shape.Jdrag_rByV,'%.3f')];
        ttl = title(titleline);
        xlab = xlabel([lineA,', ',lineB,', ',lineC]);
        set(ttl,'fontsize',titlefontsize,'interpreter','latex');
        set(xlab,'fontsize',xlabfontsize,'interpreter','latex');

        if kcase ~= 2
            colorbar('off');
        end

        % save plots
        saveas(gcf,['vfield_' num2str(kcase)],'epsc');
        saveas(gcf,['vfield_' num2str(kcase)],'png');

        %
        timecompute = toc;
        fprintf('->  Time cost plotting about %.0f seconds. \n', round(timecompute)); % about 5 min each

    end

end