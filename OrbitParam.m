classdef OrbitParam
    properties
        sma;    % [double] semi-major axis
        ecc;    % [double] eccentricity
        inc;    % [double] inclination
        lan;    % [double] longitude of ascending node
        aop;    % [double] argument of periapsis
    end
    
    
    methods
        % constructor
        function op = OrbitParam(varargin)
            switch ( nargin )
                case 0
                    op.sma = 1;
                    op.ecc = 0.1;
                    op.inc = 0.1;
                    op.lan = 0.1;
                    op.aop = 0.1;
                    
                case 5
                    
            end
        end
    end
    
    
end