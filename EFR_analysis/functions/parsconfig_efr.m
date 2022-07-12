function [ pars ] = parsconfig_efr( method )
%%
%*************************************************************************
% Copyright (c) 2015 Markus Pelz, Carl-von-Ossietzky-University Oldenburg
% Author: Markus Pelz, markus.pelz@uni-oldenburg.de
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
% ************************************************************************

%% Check input parameters are appropriate

p = inputParser;

checkmethod = @(method) (ischar(method)); % check if input is character
addRequired(p, 'method', checkmethod);

parse(p, method);
method = p.Results.method;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% general parameter

pars.fs = 16384;
pars.nchannels = 1;%1:32;
pars.totalchan = 64; % total number of used EEG channels in EEG setup
pars.nDraws = 200; % amount of draws from data for bootstrapping
%pars.nperDraw = 400;
pars.nDraws_noise = 1000;

%% condition specific parameter

switch lower(method) % change input to lower case if necessary
    
    case 'bootstrap' % Parameter for speechnoise
        
        % here is room for boostrap specific parameters    
        
    otherwise
        
        error(['There is no method called: ',method])
        
end
pars = orderfields(pars); % order fields befor outputting parameter structure

end


