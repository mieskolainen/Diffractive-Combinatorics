% Minimum and Maximum absolute variations between matrices
%
% EXAMPLE:
% 
% maxvar = zeros(size(CS{1}));
% minvar = zeros(size(CS{1}));
%
% for ID1 = 1:length(runID)
%     for ID2 = 1:length(runID)
%         [minvar, maxvar] = minmaxvar(CS{ID1}, CS{ID2}, minvar, maxvar);
%     end
% end
%
% mikael.mieskolainen@cern.ch, 2019

function [minvar, maxvar] = minmaxvar(X, Y, minvar, maxvar)

for i = 1:size(X,1)
    for j = 1:size(X,2)
        
        % Min
        value = X(i,j) - Y(i,j);
        if (value < minvar(i,j))
            minvar(i,j) = value;
        end
        
        % Max
        value = Y(i,j) - X(i,j);
        if (value > maxvar(i,j))
            maxvar(i,j) = value;
        end
    end
end

end
