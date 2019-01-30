% Create combination legends
%
% <000>
% <100>
% <010>
% ...
%
% mikael.mieskolainen@cern.ch, 2018
% Licensed under the MIT License <http://opensource.org/licenses/MIT>.

function [legends,mark] = makelegend(N, reverse)

if (nargin == 1)
    reverse = false;
end

CM = createCBM(N,reverse); % Create binary matrix

legends = cell(2^N,1);
for i = 1:2^N
    str = '';
    for k = 1:N, str = [str sprintf('%d',CM(i,k))]; end
    legends{i} = sprintf('$\\langle %s \\rangle$', str);
end
% Create linestyles
mark = {'k-'};
for i = 1:2^N-1
    mark{end+1} = '-';
end

end