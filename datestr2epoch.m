function epoch = datestr2epoch(str,F)

if nargin == 1; F=[]; end

epoch  = ( datenum(str,F) - datenum('1/1/1970') ) * (60*60*24);
