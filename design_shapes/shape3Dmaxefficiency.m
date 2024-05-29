classdef shape3Dmaxefficiency < shape3Dadjoint
    % (c) 2023 Ruowen Liu
    % This is a subclass determining optimal slip and max efficiency

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

        function obj = shape3Dmaxefficiency(design_vec)
            obj = obj@shape3Dadjoint(design_vec);
            F = zeros(obj.Nu, 1);
            trac = zeros(length(obj.t), obj.Nu);
            for k = 1:obj.Nu
                us_k = obj.W_matrix(:,k).*obj.tang;
                rhs_full_k = [real(us_k);imag(us_k)];
                mu_k = obj.A\rhs_full_k;  % solve for density
                f_k = (eye(size(obj.T))/2 + obj.T)*mu_k;   % traction
                trac(:,k) = f_k(1:end/2) + 1i*f_k(end/2+1:end);
                F(k) = 2*pi*obj.ws'*(imag(trac(:,k)).*real(obj.x));
            end

            Amatrix = zeros(obj.Nu);
            Qmatrix = zeros(obj.Nu);
            for p = 1:obj.Nu
                for q = 1:obj.Nu
                    Amatrix(p,q) = 2*pi*obj.ws'*(dotv(trac(:,p),(obj.W_matrix(:,q).*obj.tang)).*real(obj.x));
                    Qmatrix(p,q) = Amatrix(p,q) - F(p)*F(q)/obj.F0;
                end
            end

            obj.optimal_xi_slip = - Qmatrix\F; % negative sign to make U>0
            obj.uslip =  obj.W_matrix * obj.optimal_xi_slip;
            obj.JE = - F'*obj.optimal_xi_slip/obj.F0;

            [obj.trac, obj.pres, obj.U, obj.JW, obj.JD] = forwardproblem(obj);

        end % function shape3Dmaxefficiency

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
        end % plotuslipU1

        function printresults(obj)
            fprintf('-- Results:ReducedVolume %.6f,MaxEfficiency %.6f,DragForce %.6f --\n',obj.rvol,obj.JE,obj.Jdrag_rByV);
        end % printresults

        function plotblack(obj)
            figure, hold on, colorchoice = '#000000'; % black
            xzfullround = [obj.endpt0; obj.x; obj.endpt1; flipud(-real(obj.x)+1i*imag(obj.x)); obj.endpt0; obj.x(1:10)];
            plot(real(xzfullround),imag(xzfullround),'-','color',colorchoice,'linewidth',1.5);
            titlename = sprintf('$\\nu$ = %.5f, Max Efficency = %.5f, Drag Force = %.5f',obj.rvol,obj.JE,obj.Jdrag_rByV);
            tlt = title(titlename,'FontSize',15);
            set(tlt,'Interpreter','latex');
            set(gca,'FontSize',15,'TickLabelInterpreter','latex');
            grid off; axis equal; ylim([-2.5,2.5]); drawnow;
        end % plotblack

        function plotorange(obj)
            figure, hold on, colorchoice = '#F28522'; % orange
            xzfullround = [obj.endpt0; obj.x; obj.endpt1; flipud(-real(obj.x)+1i*imag(obj.x)); obj.endpt0; obj.x(1:10)];
            plot(real(xzfullround),imag(xzfullround),'-','color',colorchoice,'linewidth',1.5);
            titlename = sprintf('$\\nu$ = %.5f, Max Efficency = %.5f, Drag Force = %.5f',obj.rvol,obj.JE,obj.Jdrag_rByV);
            tlt = title(titlename,'FontSize',15);
            set(tlt,'Interpreter','latex');
            set(gca,'FontSize',15,'TickLabelInterpreter','latex');
            grid off; axis equal; ylim([-2.5,2.5]); drawnow;
        end % plotorange

        function plotgreen(obj)
            figure, hold on, colorchoice = '#77AC30'; % green
            xzfullround = [obj.endpt0; obj.x; obj.endpt1; flipud(-real(obj.x)+1i*imag(obj.x)); obj.endpt0; obj.x(1:10)];
            plot(real(xzfullround),imag(xzfullround),'-','color',colorchoice,'linewidth',1.5);
            titlename = sprintf('$\\nu$ = %.5f, Max Efficency = %.5f, Drag Force = %.5f',obj.rvol,obj.JE,obj.Jdrag_rByV);
            tlt = title(titlename,'FontSize',15);
            set(tlt,'Interpreter','latex');
            set(gca,'FontSize',15,'TickLabelInterpreter','latex');
            grid off; axis equal; ylim([-2.5,2.5]); drawnow;
        end % plotgreen

        function plotblue(obj)
            figure, hold on, colorchoice = '#0000a7'; % blue
            xzfullround = [obj.endpt0; obj.x; obj.endpt1; flipud(-real(obj.x)+1i*imag(obj.x)); obj.endpt0; obj.x(1:10)];
            plot(real(xzfullround),imag(xzfullround),'-','color',colorchoice,'linewidth',1.5);
            titlename = sprintf('$\\nu$ = %.5f, Max Efficency = %.5f, Drag Force = %.5f',obj.rvol,obj.JE,obj.Jdrag_rByV);
            tlt = title(titlename,'FontSize',15);
            set(tlt,'Interpreter','latex');
            set(gca,'FontSize',15,'TickLabelInterpreter','latex');
            grid off; axis equal; ylim([-2.5,2.5]); drawnow;
        end % plotblue

    end % methods

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

end % function forwardproblem
