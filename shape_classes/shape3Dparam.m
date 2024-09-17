classdef shape3Dparam
    % (c) 2023 Ruowen Liu
    % set initial parameters

    properties
        % NOTE:
        % p: number of points in one panel (fixed to 10)
        % np: number of panels on generating curve
        % NL: number of equi-spaced subintervals (for shape) on [0,L]
        % NLuslip: MUTST BE EVEN, number of equi-spaced subintervals (for uslip) on [0,Luslip]
        p {mustBeNumeric} = 10
        np {mustBeNumeric} = 60
        NL {mustBeNumeric} = 21 
        NLuslip {mustBeNumeric} = 200
        L {mustBeNumeric} = pi
        Luslip {mustBeNumeric} = 2*pi        
        Nu {mustBeNumeric}
    end % properties

    methods
        function obj = shape3Dparam
            obj.Nu = (obj.NLuslip-2)/2;
        end % shape3Dparam
    end % methods

end
