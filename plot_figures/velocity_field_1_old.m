% Copyright @ Ruowen Liu, July 2023

% close all
clear
warningid = 'MATLAB:nearlySingularMatrix'; warning('off',warningid);

load('/Users/ruowen/Documents/MATLAB/shapeoptswimmer2023/swimmer-con/ArbitraryInitial/result-70-2023M07D17h23m51-big-3-320/data_final.mat');

iter_total = size(design_vec_iteration,2)-1;
for iter_num = 12 % [0, 2, 4, 8, 12, iter_total]

    tic

    design_vec = design_vec_iteration(:,iter_num+1);
    shape = shape3Dmaxefficiency2(design_vec);

    % verify U is one
    disp(['U is ' num2str(shape.U,'%.4f')]);

    timeshape = toc;
    fprintf('\nCalculate shape time cost %.0f seconds. \n', round(timeshape));

    shape.printresults;

    %% prepare plot 3D swimmer view
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

    %% Calculate velocity flow
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
    rhs = [shape.uslip; shape.uslip].*[real(shape.tang);imag(shape.tang)]+[zeros(size(shape.t));shape.U*ones(size(shape.t))];   % rhs
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

    save(['iter_' num2str(iter_num) '.mat'])
    timecompute = toc;
    fprintf('->  Generating data time cost about %.0f minutes. \n', round(timecompute/60)); % about 5 min each

    %% plotting
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
    bm = brewermap([],"-Spectral");
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

    titleline = sprintf('Iteration %i',iter_num);
    lineA = ['$\nu$ ' num2str(shape.rvol,'%.3f')];lineB = ['$E$ ' num2str(shape.JE,'%.3f')];lineC = ['$J_{\! drag}$ ' num2str(shape.Jdrag_rByV,'%.3f')];
    ttl = title(titleline);
    xlab = xlabel([lineA,', ',lineB,', ',lineC]);
    set(ttl,'fontsize',titlefontsize,'interpreter','latex');
    set(xlab,'fontsize',xlabfontsize,'interpreter','latex');

    % save plots
    saveas(gcf,['iter_' num2str(iter_num)],'png');

    timecompute = toc;
    fprintf('->  Time cost plotting about %.0f minutes. \n', round(timecompute/60)); %  about 7 min for n=512, 13 min for n=700

end
