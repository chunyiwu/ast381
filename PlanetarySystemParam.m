classdef PlanetarySystemParam
%%
% PlanetarySystemParam
%
% An object holding information about a planetary system.
%
% REFERENCES
% [1] Murray, C. D., Correia, A. C. M., "Keplerian Orbits and Dynamics of 
%       Exoplanets"
% [2] Vallado, D. A., "Fundamentals of Astrodynamics and Applications"
% [3] Wikipedia, "Julian Day"
%
% AUTHOR
% Chun-Yi Wu


    properties
        sma;    % [double] semi-major axis
        ecc;    % [double] eccentricity
        inc;    % [double] inclination
        lan;    % [double] longitude of ascending node
        aop;    % [double] argument of periapsis
        tp;     % [double] time of periapsis passage (MATLAB datenum)
        
        m1;     % [double] primary body mass
        m2;     % [double] secondary body mass
        
        R1;     % [double] primary body radius
        R2;     % [double] secondary body radius
        
        r;      % [double] distance from Earth
        
        name1;  % [char] name of primary body
        name2;  % [char] name of secondary body
    end
    
    properties ( SetAccess = private )
        RM_pqw2eci;     % [double matrix 3x3] 
                        % rotation matrix from perifocal to inertial frame
    end
    
    properties ( Dependent )
        mass_ratio;     % [double] mass ratio
        mu;             % [double] gravitational parameter (=G*(m1+m2))
        n;              % [double] mean motion (=sqrt(mu/sma^3))
    end
    
    properties ( Constant )
        G = 6.67408e-11;        % gravitational constant [m^3/kg-s^2]
        
        M_Sun = 1.99855e30;     % solar mass [kg]
        M_Jup = 1.898e27;       % Jupiter mass [kg]
        M_Nep = 1.0243e26;      % Neptune mass [kg]
        M_Earth = 5.9722e24;    % Earth mass [kg]
        
        R_Sun = 6.957e8;        % solar radius [m]
        R_Jup = 69.911e6;       % Jupiter radius [m]
        R_Nep = 24.622e6;       % Neptune radius [m]
        R_Earth = 6371000;      % Earth radius [m]
        
        AU = 149597870700;      % astronomical unit [m]
        pc = 3.0857e16;         % parsec [m]
        deg = pi/180;           % radians to degrees converter
        day = 86400;            % day [s]
        yr = 86400 * 365.2422;  % year [s]
        J2000 = 2451545.0;      % Julian date at J2000
    end
    
    
    methods
        % constructor
        function op = OrbitParam(varargin)
            switch ( nargin )
                case 0
                    op.sma = OrbitParam.AU;
                    op.ecc = 0.0;
                    op.inc = 0.0;
                    op.lan = 0.0;
                    op.aop = 0.0;
                    op.tp  = OrbitParam.juldate2datenum(2451545.0);
                    op.m1  = OrbitParam.M_Sun;
                    op.m2  = OrbitParam.M_Earth;
                    op.R1  = OrbitParam.R_Sun;
                    op.R2  = OrbitParam.R_Earth;
                    op.r   = OrbitParam.pc;
                    
                    op.name1 = 'Sun';
                    op.name2 = 'Earth';
                    
                case 13
                    op.sma = varargin{ 1};
                    op.ecc = varargin{ 2};
                    op.inc = varargin{ 3};
                    op.lan = varargin{ 4};
                    op.aop = varargin{ 5};
                    op.tp  = varargin{ 6};
                    op.m1  = varargin{ 7};
                    op.m2  = varargin{ 8};
                    op.R1  = varargin{ 9};
                    op.R2  = varargin{10};
                    op.r   = varargin(11);
                    
                    op.name1 = varargin{12};
                    op.name2 = varargin{13};
                    
                otherwise
                    error('input format unrecognized');
                    
            end
            
            sO = sin(op.lan); cO = cos(op.lan);
            si = sin(op.inc); ci = cos(op.inc);
            sw = sin(op.aop); cw = cos(op.aop);
            
            RO = [ cO, sO,  0;
                  -sO, cO,  0; 
                    0,  0,  1];
            Ri = [  1,  0,  0;
                    0, ci, si;
                    0,-si, ci];
            Rw = [ cw, sw,  0;
                  -sw, cw,  0; 
                    0,  0,  1];
                
            op.RM_pqw2eci = RO'*Ri'*Rw';
        
        end
        
        
        
        % dependent properties
        function mass_ratio = get.mass_ratio(op)
            mass_ratio = op.m2 / ( op.m1 + op.m2 );
        end
        function mu = get.mu(op)
            mu = op.G * ( op.m1 + op.m2 );
        end
        function n = get.n(op)
            n = sqrt(op.mu/op.sma^3);
        end
        
        
        
        % relative state calculations
        function [r,v] = oe2rv(op,f)
        %%
        % oe2rv
        %
        % Calculate the relative position and velocity vectors from the orbital
        % elements and true anomaly.
        %
        % INPUT
        % op - [OrbitParam] orbital parameter object
        % f  - [double] true anomaly
        %
        % OUTPUT
        % r  - [double array 3x1] relative position vector
        % v  - [double array 3x1] relative velocity vector
        %
        % REFERENCE
        % [1] Murray, C. D., Correia, A. C. M., "Keplerian Orbits and Dynamics
        %       of Exoplanets"
        % [2] Vallado, D. A., "Fundamentals of Astrodynamics and
        %       Applications"
        
        p = op.sma * ( 1 - op.ecc^2 ); 
        
        rmag = p / ( 1 + op.ecc * cos(f) );
        r = rmag * op.RM_pqw2eci * [cos(f); sin(f); 0];
        
        vmag = sqrt(op.mu/p);
        v = vmag * op.RM_pqw2eci * [-sin(f); op.ecc+cos(f); 0];
        
        end
        
        
        
        % astrometry
        function [] = plotRelAstrometry(op,fnum)
        %%
        % plotRelAstrometry
        %
        % Plot the relative astrometry of the system.
        %
        % INPUT
        % fnum - [int] figure number
        
        figure(fnum); clf(fnum); hold on; axis equal; grid on; axis tight;
        
        fs = linspace(0,2*pi,3601);
        rs = zeros(3,3601);
        
        for ( i = 1 : 3601 )
            [r,~] = op.oe2rv(fs(i));
            rs(:,i) = r;
        end
        
        plot(0,0,'r+');
        plot(rs(1,:)/op.AU,rs(2,:)/op.AU,'b-');
        xlabel('X [AU]'); ylabel('Y [AU]');
        title('Relative astrometry');
        legend(op.name1,op.name2);
        end
        function [] = plotAbsAstrometry(op,fnum)
        %%
        % plotAbsAstrometry
        %
        % Plot the absolute astrometry of the system.
        %
        % INPUT
        % fnum - [int] figure number
        
        figure(fnum); clf(fnum); hold on; axis equal; grid on; axis tight;
        
        fs = linspace(0,2*pi,3601);
        rs = zeros(3,3601);
        
        for ( i = 1 : 3601 )
            [r,~] = op.oe2rv(fs(i));
            rs(:,i) = r;
        end
        
        mr1 =   - op.mass_ratio;
        mr2 = 1 - op.mass_ratio;
        
        plot(mr1*rs(1,:)/op.AU,mr1*rs(2,:)/op.AU,'r-');
        plot(mr2*rs(1,:)/op.AU,mr2*rs(2,:)/op.AU,'b-');
        plot(0,0,'k+');
        xlabel('X [AU]'); ylabel('Y [AU]');
        title('Absolute astrometry');
        legend(op.name1,op.name2,'CoM');
        end
        function [] = plotRadialVelocity(op,fnum,t0,tf)
        %%
        % plotRadialVelocity
        %
        % Plot the radial velocity of primary body.
        %
        % INPUT
        % fnum - [int] figure number
        % t0   - [double] starting time in datenum
        % tf   - [double] ending time in datenum
            
        figure(fnum); clf(fnum); hold on; grid on; axis tight;
        
        ts = linspace(t0,tf,10001);
        vs = zeros(1,10001);
        
        for ( i = 1 : 10001 )
            t = ts(i);
            
            M = mod((t-op.tp)*86400*op.n,2*pi);
            
            E = M;
            
            dE = 1000;
            while ( abs(dE) > 1e-8 )
                dE = (M-(E-op.ecc*sin(E))) / (1-op.ecc*cos(E));
                E = E + dE;
            end
            
            denom = 1 / ( 1 - op.ecc*cos(E) );
            cosf = ( cos(E) - op.ecc ) * denom;
            sinf = sqrt(1-op.ecc^2)*sin(E) * denom;
            
            f = atan2(sinf,cosf);
            
            [~,v] = op.oe2rv(f);
            vs(i) = v(3);
        end
        
        if ( max(vs)-min(vs) < 500 )
            plot(ts,vs,'b-');
            datetick('x');
            xlabel('time');
            ylabel('RV [m/s]');
            title(sprintf('Radial Velocity - %s',op.name1));
        else
            plot(ts,vs/10000,'b-');
            datetick('x');
            xlabel('time');
            ylabel('RV [km/s]');
            title(sprintf('Radial Velocity - %s',op.name1));
        end
%         datacursormode(@xaxis_time_datatip);
        end
        function [] = plotSeparation(op,fnum,t0,tf)
        %% 
        % plotSeparation
        %
        % Plot the distance between primary and secondary body during given
        % time period.
        %
        % INPUT
        % fnum - [int] figure number
        % t0   - [double] starting time in datenum
        % tf   - [double] ending time in datenum
        
        figure(fnum); clf(fnum); hold on; axis tight;
        
        ts = linspace(t0,tf,10001);
        rs = zeros(1,10001);
        zs = zeros(1,10001);
        
        for ( i = 1 : 10001 )
            t = ts(i);
            
            M = mod((t-op.tp)*86400*op.n,2*pi);
            
            E = M;
            
            dE = 1000;
            while ( abs(dE) > 1e-8 )
                dE = (M-(E-op.ecc*sin(E))) / (1-op.ecc*cos(E));
                E = E + dE;
            end
            
            denom = 1 / ( 1 - op.ecc*cos(E) );
            cosf = ( cos(E) - op.ecc ) * denom;
            sinf = sqrt(1-op.ecc^2)*sin(E) * denom;
            
            f = atan2(sinf,cosf);
            
            [r,~] = op.oe2rv(f);
            rs(i) = norm(r(1:2));
            zs(i) = r(3);
        end
        
        plot(ts,rs/op.AU,'b-');
        plot(ts,zs/op.AU,'r-');
        plot(ts,ones(size(ts))*(op.R1+op.R2)/op.AU,'k--');
        datetick('x');
        xlabel('time');
        ylabel('distance [AU]');
        legend('r=(x^2+y^2)^{1/2}','z');
        title(sprintf('Distance - %s',op.name1));
        end
    end
    
    
    
    methods ( Static )
        function JD = datenum2juldate(dnum)
        %%
        % datenum2juldate
        %
        % Calculate the Julian date from the MATLAB datenum.
        %
        % INPUT
        % dnum - [double] MATLAN datenum
        %
        % OUTPUT
        % JD - [double] Julian day
        %
        % REFERENCE
        % [1] Wikipedia, "Julian Day"
        
        dv = datevec(dnum);
        
        year = dv(1);
        month = dv(2);
        day = dv(3);
        hour = dv(4);
        minute = dv(5);
        second = dv(6);
        
        a = floor((14-month)/12);
        y = year + 4800 - a;
        m = month + 12 * a - 3;
        
        JDN = day + floor( (153*m+2)/5 ) + 365*y + floor(y/4) - floor(y/100) + floor(y/400) - 32045;
        JD = JDN + (hour-12)/24 + minute/1440 + second/86400;
        end
        function dnum = juldate2datenum(JD)
        %%
        % juldate2datenum
        %
        % Calcualte datenum from Julian day.
        %
        % INPUT
        % JD - [double] Julian day
        %
        % OUTPUT
        % dnum - [double] MATLAB datenum
        %
        % REFERENCE
        % [1] Wikipedia, "Julian Day"
        
        % constants
        y = 4716; j = 1401; m = 2; n = 12; r = 4; p = 1461;
        v = 3; u = 5; s = 153; w = 2; B = 274277; C = -38;
        
        % fractional day
        t_frac = (JD - round(JD))*86400 + 43200;
        S = mod(t_frac,60); t_frac = (t_frac-S)/60;
        M = mod(t_frac,60); 
        H = floor(t_frac/60);
        
        JD = round(JD);
        f = JD + j + floor(((floor((4*JD+B)/146097))*3)/4)+C;
        e = r*f + v;
        g = floor(mod(e,p)/r);
        h = u*g + w;
        D = floor((mod(h,s))/u) + 1;
        Mo = mod(floor(h/s)+m,n) + 1;
        Y = floor(e/p) - y + floor((n+m-Mo)/n);
        
        dnum = datenum(Y,Mo,D,H,M,S);
        end
    end
end