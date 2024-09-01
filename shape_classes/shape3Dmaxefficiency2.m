classdef shape3Dmaxefficiency2 < shape3Dadjoint
    % This is a subclass determining optimal slip and max efficiency
    % use new slip optimization method
    % (c) 2023 Ruowen Liu

    properties
        uslip {mustBeNumeric}
        optimal_xi_slip {mustBeNumeric}
        JE {mustBeNumeric}
        trac {mustBeNumeric}
        pres {mustBeNumeric}
        U {mustBeNumeric}
        JW {mustBeNumeric}
        JD {mustBeNumeric}
    end % properties

    methods

        function obj = shape3Dmaxefficiency2(design_vec)
            obj = obj@shape3Dadjoint(design_vec); % adjoint solution
            Mat = [ (real(obj.tang).*obj.C(1:end/2,:) + imag(obj.tang).*obj.C(end/2+1:end,:)) ; ...
                (real(obj.nx).*obj.A(1:end/2,:) + imag(obj.nx).*obj.A(end/2+1:end,:)) ];
            rhs = [dotv(obj.trachat, obj.tang) ; zeros(size(obj.nx))];
            mu0 = Mat\rhs;
            trac_til = obj.C(1:end/2,:)*mu0 + 1i*obj.C(end/2+1:end,:)*mu0;
            u_til = obj.A(1:end/2,:)*mu0 + 1i*obj.A(end/2+1:end,:)*mu0;
            zs = dotv(u_til, obj.tang);
            Mat_full = [obj.A [zeros(length(obj.t),1);-ones(length(obj.t),1)];...
                (real(obj.x).*obj.ws(:))'*obj.C(end/2+1:end,:) 0];
            rhs_full = [real(u_til);imag(u_til);0];
            sol_full = Mat_full\rhs_full;
            mu0 = sol_full(1:end-1);
            U = sol_full(end);  % solve for U
            obj.JE = -U/(1+U);
            %normalize U
            obj.uslip = zs/U;
            [obj.trac, obj.pres, obj.U, obj.JW, obj.JD] = forwardproblem(obj);
        end % shapeDerivative

        function plotuslip(obj)
            figure
            plot(obj.arclen/obj.ell, obj.uslip);
            xlabel('$s/\ell$','Interpreter','latex');
            ylabel('$u^{{S}}$','Interpreter','latex');
            set(gca, 'fontsize', 20);
            xlim([0,1]);
            ylim([0,6]);
            grid on;
            drawnow;
        end % plotuslip

        function plotuslipU1(obj)
            figure
            plot(obj.arclen/obj.ell, obj.uslip./obj.U, 'LineWidth', 1.5);
            xlabel('$s/\ell$','Interpreter','latex');
            ylabel('$u^{{S}}$ (reset $U=1$)','Interpreter','latex');
            set(gca, 'fontsize', 20);
            xlim([0,1]);
            ylim([0,1.5]);
            grid on;
            drawnow;
        end % plotuslip

        function printresults(obj)
            fprintf('-- Results:ReducedVolume %.6f,MaxEfficiency %.6f,DragForce %.6f --\n',obj.rvol,obj.JE,obj.Jdrag_rByV);
        end % printresults

        function plotblack(obj)
            figure, hold on, colorchoice = '#000000'; % black
            xzfullround = [obj.endpt0; obj.x; obj.endpt1; flipud(-real(obj.x)+1i*imag(obj.x)); obj.endpt0; obj.x(1:10)];
            plot(real(xzfullround),imag(xzfullround),'-','color',colorchoice,'linewidth',1.5);
            titlename = sprintf('$\\nu$ = %.3f, $E$ = %.3f, $J_{drag}$ = %.3f',obj.rvol,obj.JE,obj.Jdrag_rByV);
            tlt = title(titlename,'FontSize',15);
            set(tlt,'Interpreter','latex');
            set(gca,'FontSize',15,'TickLabelInterpreter','latex');
            grid off; axis equal; ylim([-2.5,2.5]); drawnow;
        end

        function plotorange(obj)
            figure, hold on, colorchoice = '#F28522'; % orange
            xzfullround = [obj.endpt0; obj.x; obj.endpt1; flipud(-real(obj.x)+1i*imag(obj.x)); obj.endpt0; obj.x(1:10)];
            plot(real(xzfullround),imag(xzfullround),'-','color',colorchoice,'linewidth',1.5);
            titlename = sprintf('$\\nu$ = %.5f, Max Efficency = %.5f, Drag Force = %.5f',obj.rvol,obj.JE,obj.Jdrag_rByV);
            tlt = title(titlename,'FontSize',15);
            set(tlt,'Interpreter','latex');
            set(gca,'FontSize',15,'TickLabelInterpreter','latex');
            grid off; axis equal; ylim([-2.5,2.5]); drawnow;
        end

        function plotgreen(obj)
            figure, hold on, colorchoice = '#77AC30'; % green
            xzfullround = [obj.endpt0; obj.x; obj.endpt1; flipud(-real(obj.x)+1i*imag(obj.x)); obj.endpt0; obj.x(1:10)];
            plot(real(xzfullround),imag(xzfullround),'-','color',colorchoice,'linewidth',1.5);
            titlename = sprintf('$\\nu$ = %.5f, Max Efficency = %.5f, Drag Force = %.5f',obj.rvol,obj.JE,obj.Jdrag_rByV);
            tlt = title(titlename,'FontSize',15);
            set(tlt,'Interpreter','latex');
            set(gca,'FontSize',15,'TickLabelInterpreter','latex');
            grid off; axis equal; ylim([-2.5,2.5]); drawnow;
        end

        function plotblue(obj)
            figure, hold on, colorchoice = '#0000a7'; % blue
            xzfullround = [obj.endpt0; obj.x; obj.endpt1; flipud(-real(obj.x)+1i*imag(obj.x)); obj.endpt0; obj.x(1:10)];
            plot(real(xzfullround),imag(xzfullround),'-','color',colorchoice,'linewidth',1.5);
            titlename = sprintf('$\\nu$ = %.5f, Max Efficency = %.5f, Drag Force = %.5f',obj.rvol,obj.JE,obj.Jdrag_rByV);
            tlt = title(titlename,'FontSize',15);
            set(tlt,'Interpreter','latex');
            set(gca,'FontSize',15,'TickLabelInterpreter','latex');
            grid off; axis equal; ylim([-2.5,2.5]); drawnow;
        end

    end

end

function [trac, pres, U, JW, JD] = forwardproblem(s)

Mat_full = [s.A [zeros(length(s.t),1);-ones(length(s.t),1)];...
    (real(s.x).*s.ws(:))'*s.C(end/2+1:end,:) 0];

us = s.uslip.*s.tang; % only (optimal) slip (U excluded)

rhs_full = [real(us);imag(us);0];
sol_full = Mat_full\rhs_full;
mu0 = sol_full(1:end-1);
U = sol_full(end);  % solve for U
traction = s.C*mu0;

trac = traction(1:end/2) + 1i*traction(end/2+1:end);
pres = s.P*mu0-1/2*(real(s.nx).*mu0(1:end/2)+imag(s.nx).*mu0(end/2+1:end));
JW = 2*pi*(real(s.x).*s.ws(:))'*dotv(trac, s.uslip.*s.tang+1i*U);
JD = s.F0*U*U;
% JE = JD/JW;

end % function
