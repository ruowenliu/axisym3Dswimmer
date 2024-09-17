classdef shape3Dbasic_resetNL
    % (c) 2024 Ruowen Liu
    % This is an object-oriented programming for basic shape properties
    % Use it as the superclass for all further computations

    properties
        % NOTE:
        % p: number of points in one panel (fixed to 10)
        % np: number of panels on generating curve
        % NL: number of equi-spaced subintervals (for shape) on [0,L]
        % NLuslip: MUTST BE EVEN, number of equi-spaced subintervals (for uslip) on [0,Luslip]
        p {mustBeNumeric} = 10
        np {mustBeNumeric} = 60
        NL {mustBeNumeric}
        NLuslip {mustBeNumeric} = 200
        L {mustBeNumeric} = pi
        Luslip {mustBeNumeric} = 2*pi        
        Nu {mustBeNumeric}

        xi_R {mustBeNumeric}
        xi_Z {mustBeNumeric}
        design_vec {mustBeNumeric}
        tpan {mustBeNumeric}
        tlo {mustBeNumeric}
        thi {mustBeNumeric}
        xlo {mustBeNumeric}
        xhi {mustBeNumeric}
        w {mustBeNumeric}
        ws {mustBeNumeric}
        x {mustBeNumeric}
        xp {mustBeNumeric}
        xpp {mustBeNumeric}
        sp {mustBeNumeric}
        tang {mustBeNumeric}
        nx {mustBeNumeric}
        cur {mustBeNumeric}
        t {mustBeNumeric}
        wxp {mustBeNumeric}
        arclen {mustBeNumeric}
        ell {mustBeNumeric}
        vol {mustBeNumeric}
        area {mustBeNumeric}
        rvol {mustBeNumeric}
        endpt0 {mustBeNumeric}
        endpt1 {mustBeNumeric}
        W_matrix {mustBeNumeric}
        tL {mustBeNumeric}
        XL {mustBeNumeric}
    end % properties

    methods

        function obj = shape3Dbasic_resetNL(NL, design_vec)
            obj.NL = NL;
            obj.Nu = (obj.NLuslip-2)/2;
            obj.design_vec = design_vec;
            obj.xi_R = [0;design_vec(1:end/2);0];
            obj.xi_Z = [0;design_vec((1+(end/2)):end);0];
            % % reset the first and last coefficients to make pole nice
            % % x(0)=x(pi)=0
            % fprintf('old xi_R(1) xi_R(end) %f %f \n',obj.xi_R(1),obj.xi_R(end));
            obj.xi_R(1) = -[26,66,26,1]*obj.xi_R(2:5);
            obj.xi_R(end) = -[1,26,66,26]*obj.xi_R(end-4:end-1);
            % fprintf('new xi_R(1) xi_R(end) %f %f \n',obj.xi_R(1),obj.xi_R(end));
            % % z'(0)=z'(pi)=0
            % fprintf('old xi_Z(1) xi_Z(end) %f %f \n',obj.xi_Z(1),obj.xi_Z(end));
            obj.xi_Z(1) = ([-10,0,10,1]*obj.xi_Z(2:5));
            obj.xi_Z(end) = -([-1,-10,0,10]*obj.xi_Z(end-4:end-1));
            % fprintf('new xi_Z(1) xi_Z(end) %f %f \n',obj.xi_Z(1),obj.xi_Z(end));

            obj.tpan = (0:(obj.L/obj.np):obj.L)'; % for panels
            str.Z = @(t) B0_original( (t*obj.NL/obj.L) - ((-5):(obj.NL-1)) ) * ( obj.xi_R + 1i * obj.xi_Z );
            str.p = obj.p; str.tpan = obj.tpan; [str,~] = quadrp(str, obj.p*obj.np, 'p', 'G');
            obj.t = str.t; obj.tlo = str.tlo; obj.thi = str.thi; obj.xlo = str.xlo; obj.xhi = str.xhi;
            obj.x = str.x; obj.xp = str.xp; obj.xpp = str.xpp; obj.sp = str.sp; obj.nx = str.nx; obj.tang = str.tang;
            obj.w = str.w; obj.wxp = str.wxp; obj.ws = str.ws;
            obj.cur = str.cur; obj.ell = sum(obj.ws);

            % location of BSP control points
            obj.tL = linspace(0,pi,obj.NL+1)';
            obj.XL = str.Z(obj.tL);

            % % add a negative sign because the curve is going downward, so dy/dt flips
            % % Vy = int_{a}^{b} pi x^2 (dy/dt) dt
            % % https://en.wikipedia.org/wiki/Solid_of_revolution
            obj.vol = -pi*obj.w'*((real(obj.x).^2).*imag(obj.xp));

            % % surfacearea = 2pi \int_0^{ell} R ds
            obj.area = 2*pi*obj.ws'*real(obj.x);

            % % reduced_volume = 6sqrt(pi)V/(A^(2/3))
            obj.rvol = 6*sqrt(pi)*obj.vol/(obj.area^1.5);

            obj.arclen = zeros(size(obj.t));
            len_obj = length(obj.t);
            for k = 1:len_obj % cannot use parfor
                tEnd = obj.t(k); idx = (obj.thi<tEnd); arcLen = sum(obj.ws(1:obj.p*sum(idx))); % entire panel \subset [0,t];
                tlo_arclen = obj.tlo(sum(idx)+1); thi_arclen = tEnd; [x0, w0, D0] = gauss(obj.p);
                t_arclen = tlo_arclen+(1+x0)/2*(thi_arclen-tlo_arclen); w_arclen = w0*(thi_arclen-tlo_arclen)/2;
                x_arclen = str.Z(t_arclen); xp_arclen = D0*x_arclen*2/(thi_arclen-tlo_arclen);
                ws_arclen = w_arclen.*abs(xp_arclen);
                obj.arclen(k) = arcLen+sum(ws_arclen);
            end

            obj.endpt0 = str.Z(0); obj.endpt1 = str.Z(obj.L);

            % W_matrix is dependent of shape. W(arclen/ell) for basis
            obj.W_matrix = find_W_matrix(obj);

        end % function shape3Dbasic

        function visualize(obj)
            figure
            hold on
            % % blue '#0072BD' % red '#A2142F' % green '#77AC30'
            plot(real([obj.endpt0;obj.x;obj.endpt1]),imag([obj.endpt0;obj.x;obj.endpt1]),'-o','color','#77AC30','MarkerSize',2,'linewidth',1.5);
            plot(-real(obj.x),imag(obj.x),'-','color','cyan','linewidth',1.5);
            plot(real(obj.x(1)),imag(obj.x(1)),'or'), text(real(obj.x(1)),imag(obj.x(1))-0.1,'start');
            grid_jump = 15; grid_span = [1:grid_jump:length(obj.x)-1, length(obj.x)];
            quiver(real(obj.x(grid_span)),imag(obj.x(grid_span)),real(obj.nx(grid_span)),imag(obj.nx(grid_span)),'color','#A2142F');
            quiver(real(obj.x(grid_span)),imag(obj.x(grid_span)),real(obj.tang(grid_span)),imag(obj.tang(grid_span)),'color','#0072BD');
            plot(real(obj.XL), imag(obj.XL), '.r', 'markersize', 10)
            set(gca, 'fontsize', 20);
            axis equal;
            ylim([-2.5,2.5]);
            grid on;
            drawnow;
        end % visualize

        function visulaize_simple(obj)
            figure, hold on, colorchoice = '#000000'; % black
            plot(real([obj.endpt0;obj.x;obj.endpt1]),imag([obj.endpt0;obj.x;obj.endpt1]),'-','color',colorchoice,'linewidth',1.5);
            plot(-real([obj.endpt0;obj.x;obj.endpt1]),imag([obj.endpt0;obj.x;obj.endpt1]),'-','color',colorchoice,'linewidth',1.5);
            titlename = sprintf('$\\nu$ = %.5f',obj.rvol);
            tlt = title(titlename,'FontSize',15);
            set(tlt,'Interpreter','latex');
            set(gca,'FontSize',15,'TickLabelInterpreter','latex');
            grid off; axis equal; ylim([-2.5,2.5]); drawnow;
        end % visualize_simple

    end % methods

end

function W = find_W_matrix(s)
% This computes the matrix W such that
% [w(tau)_1, ..., w(tau)_{NLuslip+5}] = W * xi_uslip
% Here tau = s/l is between 0 and 1

v_bas = [zeros(1,s.Nu); eye(s.Nu); zeros(1,s.Nu); -flipud(eye(s.Nu)); zeros(1,s.Nu)];

rhs = [zeros(4,s.Nu); v_bas];

MAT = zeros(size(rhs,1));
MAT(1,1:5) = [-1,-10,0,10,1]/24; MAT(1,end-4:end) = -MAT(1,1:5);
MAT(2,1:5) = [1,2,-1,2,1]/6; MAT(2,end-4:end) = -MAT(2,1:5);
MAT(3,1:5) = [-0.5,1,0,-1,0.5]; MAT(3,end-4:end) = -MAT(3,1:5);
MAT(4,1:5) = [1,-4,6,-4,1]; MAT(4,end-4:end) = -MAT(4,1:5);
for k=5:(size(rhs,1)), MAT(k,k-4:k) = [1,26,66,26,1]/120; end

c_bas = MAT \ rhs;

W = B0_original( ((s.arclen/s.ell)*s.NLuslip/2) - ((-5):(s.NLuslip-1)) ) * c_bas;

end % find_W_matrix
