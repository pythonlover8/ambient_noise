function date_str = epoch2datestr(etime,F)

if nargin == 1; F='mm/dd/yyyy HH:MM:SS.FFF'; end

% etime = (datenum(dates) -datenum('1/1/1970')) * (60*60*24);

etime = (etime./(60*60*24)) + datenum('1/1/1970') ;

date_str = datestr(datevec(etime),F);