% Median Absolute Deviation
%
% mikael.mieskolainen@cern.ch, 2018
% Licensed under the MIT License <http://opensource.org/licenses/MIT>.

function cost = madcost(residual)

cost = median(abs(residual - median(residual) ) );

end