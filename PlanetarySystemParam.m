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
% [4] "Positional Astronomy, Annual parallax"
%    <http://star-www.st-and.ac.uk/~fv/webnotes/>
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
        
        rStar;  % [double] distance from Sun
        RA;     % [double] right ascension
        dec;    % [double] declination
        
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
        J2000 = PlanetarySystemParam.juldate2datenum(2451545.0);      
                                % datenum at J2000 epoch
    end
    
    
    methods
        % constructor
        function op = PlanetarySystemParam(varargin)
            switch ( nargin )
                case 0
                    op.sma = PlanetarySystemParam.AU;
                    op.ecc = 0.0;
                    op.inc = 0.0;
                    op.lan = 0.0;
                    op.aop = 0.0;
                    op.tp  = PlanetarySystemParam.juldate2datenum(2451545.0);
                    op.m1  = PlanetarySystemParam.M_Sun;
                    op.m2  = PlanetarySystemParam.M_Earth;
                    op.R1  = PlanetarySystemParam.R_Sun;
                    op.R2  = PlanetarySystemParam.R_Earth;
                    op.rStar  = PlanetarySystemParam.pc;
                    op.RA     = 0;
                    op.dec    = 0;
                    
                    op.name1 = 'Sun';
                    op.name2 = 'Earth';
                    
                case 15
                    op.sma      = varargin{ 1};
                    op.ecc      = varargin{ 2};
                    op.inc      = varargin{ 3};
                    op.lan      = varargin{ 4};
                    op.aop      = varargin{ 5};
                    op.tp       = varargin{ 6};
                    
                    op.m1       = varargin{ 7};
                    op.m2       = varargin{ 8};
                    
                    op.R1       = varargin{ 9};
                    op.R2       = varargin{10};
                    
                    op.rStar    = varargin{11};
                    op.RA       = varargin{12};
                    op.dec      = varargin{13};
                    
                    op.name1    = varargin{14};
                    op.name2    = varargin{15};
                    
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
        function [r,v] = oe2rv_t(op,t)
        %%
        % oe2rv
        %
        % Calculate the relative position and velocity vectors from the orbital
        % elements and true anomaly.
        %
        % INPUT
        % op - [PlanetarySystemParam] planetary system parameter object
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
        
        % solve mean anomaly
        M = mod((t-op.tp)*86400*op.n,2*pi);
           
        % solve eccentric anomaly
        E = M;

        dE = 1000;
        while ( abs(dE) > 1e-8 )
            dE = (M-(E-op.ecc*sin(E))) / (1-op.ecc*cos(E));
            E = E + dE;
        end

        % solve true anomaly
        denom = 1 / ( 1 - op.ecc*cos(E) );
        cosf = ( cos(E) - op.ecc ) * denom;
        sinf = sqrt(1-op.ecc^2)*sin(E) * denom;

        f = atan2(sinf,cosf);
        
        % solve position/velocity vector
        p = op.sma * ( 1 - op.ecc^2 ); 
        
        rmag = p / ( 1 + op.ecc * cos(f) );
        r = rmag * op.RM_pqw2eci * [cos(f); sin(f); 0];
        
        vmag = sqrt(op.mu/p);
        v = vmag * op.RM_pqw2eci * [-sin(f); op.ecc+cos(f); 0];
        
        end
        function [r,v] = oe2rv(op,f)
        %%
        % oe2rv
        %
        % Calculate the relative position and velocity vectors from the orbital
        % elements and true anomaly.
        %
        % INPUT
        % op - [PlanetarySystemParam] planetary system parameter object
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
        function [r,z] = findSeparation(op,t)
        %%
        % findSeparation
        %
        % Find the separation distance between primary and secondary
        % bodies.
        
        pos = op.oe2rv_t(t);
        r = norm(pos(1:2));
        z = pos(3);
        end
        function [rv] = findRadialVelocity(op,t)
            [~,vel] = op.oe2rv_t(t);
            rv = vel(3);
        end
        
        
        
        % transit calculation
        function t = solveTransitThreshold(op,t1,t2)
        %%
        % solveTransitThreshold
        %
        % Solve for either when the transit start or end, i.e. when the
        % distance between the center of two bodies is the same as sum of
        % body radii.
        %
        % Uses bisection method, so requires a range of time 
        % with only one crossing.
        %
        % INPUT
        % op - [PlanetarySystemParam] planetary system parameter object
        % t1 - [double] lower bound for the range
        % t2 - [double] upper bound for the range
        
        %% set up the function to solve
        r = @(t) ( op.findSeparation(t) - (op.R1 + op.R2) );
        
        %% solve when r = 0
        % check if one of the bound is already satisfying the condition
        r1 = r(t1); r2 = r(t2);
        if ( abs(r1) <= 1e-2 )
            t = t1;
            return;
        end
        
        if ( abs(r2) <= 1e-2 )
            t = t2;
            return;
        end
        
        % modified NR method
        dt = t2 - t1;
        dr = r2 - r1;
        iter = 0;
        
        while ( (dt > 1e-9 && abs(dr) > 1e-3) && iter < 100 )
            iter = iter + 1;
            
            tm = (t1+t2)/2;
            rm = r(tm);
            
            if ( r1 * rm > 0 )
                t1 = tm;
                r1 = rm;
            else
                t2 = tm;
                r2 = rm;
            end
            
            dt = t2 - t1;
            dr = r2 - r1;
        end
        
        t = (t1+t2)/2;
        end
        function t = solveFrontBackThreshold(op,t1,t2)
        %% 
        % solveFrontBackThreshold
        %
        % Solve when the planet goes from being in front to behind (or vice
        % versa) the star.
        %
        % INPUT
        % op - [PlanetarySystemParam] planetary system parameter object
        % t1 - [double] lower bound for the range
        % t2 - [double] upper bound for the range
        
        %% solve when r = 0
        % check if one of the bound is already satisfying the condition
        [~,r1] = op.findSeparation(t1);
        [~,r2] = op.findSeparation(t2);
        if ( abs(r1) <= 1e-2 )
            t = t1;
            return;
        end
        
        if ( abs(r2) <= 1e-2 )
            t = t2;
            return;
        end
        
        % modified NR method
        dt = t2 - t1;
        dr = r2 - r1;
        iter = 0;
        
        while ( dt > 1e-9 && abs(dr) > 1e-3 )
            iter = iter + 1;
            
            tm = (t1+t2)/2;
            [~,rm] = op.findSeparation(tm);
            
            if ( r1 * rm > 0 )
                t1 = tm;
                r1 = rm;
            else
                t2 = tm;
                r2 = rm;
            end
            
            dt = t2 - t1;
            dr = r2 - r1;
            
            
            fprintf('%s (%g) ~ %s (%g)\n',datestr(t1,31),r1,datestr(t2,31),r2);
        end
        
        t = (t1+t2)/2;
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
            vs(i) = op.findRadialVelocity(ts(i));
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
            [r,z] = op.findSeparation(ts(i));
            rs(i) = r;
            zs(i) = z;
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
        function [] = plotRAdecGraph(op,fnum,ts,parallax,pm)
        %%
        % plotRAdecGraph
        %
        % Plot the right ascension-declination graph of the star.
        %
        % INPUT
        % fnum - [int] figure number
        % t0   - [double] starting time in datenum
        % tf   - [double] ending time in datenum
        % parallax - [logical] if considering annual parallax motion (false
        %       if omitted)
        % propmotion - [double array] proper motion of star in milliarcsecond 
        %           per year ([RA,dec]) (0 if omitted)
        %
        % REFERENCE
        % [1] "Positional Astronomy, Annual parallax"
        %    <http://star-www.st-and.ac.uk/~fv/webnotes/>
        
        %% calculations        
        tdiv = length(ts);
        ras = zeros(size(ts));
        decs = zeros(size(ts));
        year = datevec(ts(1));
        year = year(1);
        
        for ( i = 1 : tdiv )
            % heliocentric coordinate
            a = op.RA;
            d = op.dec;
            
            % add planet effect
            pos = - op.mass_ratio * op.oe2rv_t(ts(i));
            a = a - pos(2) / op.rStar;
            d = d + pos(1) / op.rStar;
            
            % add proper motion
            if ( exist('pm','var') )
                a = a + pm(1) / 3600000 / 365.2425 * (ts(i)-op.J2000);
                d = d + pm(2) / 3600000 / 365.2425 * (ts(i)-op.J2000);
            end
            
            % add parallax
            if ( parallax )
                e = op.dms2deg(23,26,00);
                p = 1/(op.rStar/op.pc) / 3600;
                ls = - 2*pi/365.2425 * (ts(i)-vernalEquinox(year));
                
                sind = sin(d); cosd = cos(d);
                sina = sin(a); cosa = cos(a);
                sine = sin(e); cose = cos(e);
                
                sinb = sind*cose - cosd*sine*sina;
                cosb = sqrt(1-sinb^2);
                sinl = (sind-sinb*cose) / (cosb*sine);
                cosl = (cosa*cosd) / (cosb);
                l = atan2(sinl,cosl);
                
                dl =  p * sin(ls-l) / cosb;
                db = -p * cos(ls-l) * sinb;
                
                a = a + dl;
                d = d + db;                
            end
            
            ras(i) = a;
            decs(i) = d;
        end
        
        %% plot
        figure(fnum); clf(fnum); 
        
        % RA versus time
        subplot(2,2,4);
        plot((ras-op.RA)*3600000,ts,'b+');
        datetick('y'); grid on;
        xlabel('\DeltaRA [mas]');
        
        % dec verus time
        subplot(2,2,1);
        plot(ts,(decs-op.dec)*3600000,'b+');
        datetick('x'); grid on;
        ylabel('\Deltadec [mas]');
        
        % RA versus dec
        subplot(2,2,2); axis equal;
        plot((ras-op.RA)*3600000,(decs-op.dec)*3600000,'b+'); grid on;
        
        % notes
        subplot(2,2,3); axis off;
        text(0,3/4,'Reference position (J2000):','Unit','Normalized');
        [h,m,s] = op.deg2hms(op.RA);
        str = sprintf('RA: %d^h %d^m %g^s',h,m,s);
        text(0.1,2/4,str,'Unit','Normalized');
        [d,m,s] = op.deg2dms(op.dec);
        str = sprintf('dec: %d^o %d'' %g"',d,m,s);
        text(0.1,1/4,str,'Unit','Normalized');
        end
    end
    
    
    
    methods ( Static )
        %% date conversion
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
        
        
        
        %% angle conversion
        function deg = hms2deg(h,m,s)
        %%
        % hms2deg
        %
        % Convert hour angle to degrees. If an entry is omitted, it's
        % assumed to be zero.
        %
        % INPUT
        % h - [double] hour 
        % m - [double] minute
        % s - [double] second
        %
        % OUTPUT
        % deg - [double] degrees
        
        if ( ~exist('m','var') )
            m = 0;
        end
        
        if ( ~exist('s','var') )
            s = 0;
        end
        
        deg = 15 * ( h + m/60 + s/3600 );
        end
        function [h,m,s] = deg2hms(deg)
        %%
        % deg2hms
        %
        % Convert degrees to hour-minute-seconds.
        %
        % INPUT
        % deg - [double] degrees
        %
        % OUTPUT
        % h - [double] hour 
        % m - [double] minute
        % s - [double] second (this is the only one that would be
        %           non-integer)
        
        deg = deg / 15;
        h = floor(deg); deg = (deg - h) * 60;
        m = floor(deg);
        s = (deg - m) * 60;
        end
        function deg = dms2deg(d,m,s)
        %%
        % dms2deg
        %
        % Convert degree-minutes-seconds to degrees. If an entry is omitted,
        % it's assumed to be zero.
        %
        % INPUT
        % h - [double] hour 
        % m - [double] minute
        % s - [double] second
        %
        % OUTPUT
        % deg - [double] degrees
        
        if ( ~exist('m','var') )
            m = 0;
        end
        
        if ( ~exist('s','var') )
            s = 0;
        end
        
        deg = ( d + m/60 + s/3600 );
        end
        function [h,m,s] = deg2dms(deg)
        %%
        % deg2dms
        %
        % Convert degrees to degree-minute-seconds.
        %
        % INPUT
        % deg - [double] degrees
        %
        % OUTPUT
        % h - [double] hour 
        % m - [double] minute
        % s - [double] second (this is the only one that would be
        %           non-integer)
        
        h = floor(deg); deg = (deg - h) * 60;
        m = floor(deg);
        s = (deg - m) * 60;
        end
        
    end
end