function dnum = vernalEquinox(year)
%%
% vernalEquinox
%
% Return the vernal equinox time in MATLAB datenum. This is a database for
% 2016 to 2020.
%
% INPUT
% year - [int] year
%
% OUTPUT
% dnum - [double] the datenum of vernal equinox
%
% AUTHOR
% Chun-Yi Wu

switch ( year )
    case 2016
        dnum = datenum(2016,3,20,04,37,00);
    case 2017
        dnum = datenum(2017,3,20,10,26,00);
    case 2018
        dnum = datenum(2018,3,20,16,16,00);
    case 2019
        dnum = datenum(2019,3,20,22,48,00);
    case 2020
        dnum = datenum(2020,3,20,03,54,00);
    otherwise
        error('currently only support 2016-2020.');   
    
end

end