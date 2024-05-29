classdef shape3Dadjoint < shape3Dbasic
    % (c) 2023 Ruowen Liu
    % This is a subclass of shape3Dbasic
    % solves the adjoint problem

    properties
        mu {mustBeNumeric} = 1 % dynamic viscosity
        A {mustBeNumeric}
        T {mustBeNumeric}
        C {mustBeNumeric}
        P {mustBeNumeric}
        F0 {mustBeNumeric}
        trachat {mustBeNumeric}
        preshat {mustBeNumeric}
        Jdrag_rByV {mustBeNumeric}
    end % properties

    methods

        function obj = shape3Dadjoint(design_vec)
            obj = obj@shape3Dbasic(design_vec);
            obj.A = AlpertSphereSLPMat(obj);
            obj.T = AlpertSphereSLPMatT(obj);
            obj.C = (eye(size(obj.T))/2 + obj.T);
            obj.P = ChunkSphereSLPMatP(obj);
            uD_adj = [0*ones(size(obj.t));ones(size(obj.t))];   % because uD=e_Z=1i
            mu_adj = obj.A\uD_adj;  % solve for density
            f_adj = obj.C*mu_adj;   % traction
            obj.trachat = f_adj(1:end/2) + 1i*f_adj(end/2+1:end);
            obj.preshat = obj.P*mu_adj-1/2*(real(obj.nx).*mu_adj(1:end/2)+imag(obj.nx).*mu_adj(end/2+1:end));
            obj.F0 = 2*pi*obj.ws'*(imag(obj.trachat).*real(obj.x));
            obj.Jdrag_rByV = obj.F0/6/pi/obj.mu/((3*obj.vol/4/pi)^(1/3));
        end % function shape3Dadjoint

        function printresults(obj)
            fprintf('-- Results: Reduced_Vol = %.6f  Drag_rByV = %.6f --\n',... 
                obj.rvol, obj.Jdrag_rByV);
        end % printresults

    end % methods

end
