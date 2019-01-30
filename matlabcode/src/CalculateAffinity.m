% mikael.mieskolainen@cern.ch, 2018
% Licensed under the MIT License <http://opensource.org/licenses/MIT>.

function [affinity] = CalculateAffinity(data)

% set the parameters
sigma = 1;

affinity = zeros(size(data,1));

for i=1:size(data,1)    
    for j=1:size(data,1)
        dist = poincaredist(data(i,:), data(j,:));
        affinity(i,j) = dist;
        %affinity(i,j) = exp(-dist/(2*sigma^2));
    end
end


end