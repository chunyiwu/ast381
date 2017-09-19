clear; clc; close all hidden;

sma = 0.449 * PlanetarySystemParam.AU;
ecc = 0.9336;
inc = 89.285 * PlanetarySystemParam.deg;
lan = 160.98 * PlanetarySystemParam.deg;
aop = 300.651 * PlanetarySystemParam.deg;
tp  = PlanetarySystemParam.juldate2datenum(2454424.857);

P = 111.43637 * PlanetarySystemParam.day;
m1 = 4*pi^2*sma^3 / P^2 / PlanetarySystemParam.G;

m2 = 3.94 * PlanetarySystemParam.M_Jup;

R1 = 0.98 * PlanetarySystemParam.R_Sun;
R2 = 0.098 * R1;

rStar = 58.4 * PlanetarySystemParam.pc;
RA  = PlanetarySystemParam.hms2deg(09,22,37.0);
dec = PlanetarySystemParam.dms2deg(50,36,13);


op = PlanetarySystemParam(sma,ecc,inc,lan,aop,tp,m1,m2,R1,R2,rStar,RA,dec, ...
    'HD 80606','HD 80606 b');

% % op.plotRelAstrometry(1);
% % op.plotAbsAstrometry(2);
% 
% t0 = datenum(2017,08,01,00,00,00);
% tf = datenum(2017,12,31,00,00,00);
% 
% % op.plotRadialVelocity(3,t0,tf);
% op.plotSeparation(4,t0,tf);
% 
t_0 = datenum(2017,12,01,00,00,00);
t_f = datenum(2017,12,15,00,00,00);
op.plotSeparation(5,t_0,t_f);

% t1 = datenum(2017,12,07,11,00,00); 
% t2 = datenum(2017,12,07,13,00,00);
% datevec(op.solveTransitThreshold(t1,t2))
% 
% t1 = datenum(2017,12,07,22,00,00);
% t2 = datenum(2017,12,08,02,00,00);
% datevec(op.solveFrontBackThreshold(t1,t2))

fnum = 6;
ts = datenum(2017,08,01) + sort(rand(100,1)) * 153;
parallax = true;
pm = [40.64,10.17];
% pm = [0,0];
op.plotRAdecGraph(fnum,ts,parallax,pm)