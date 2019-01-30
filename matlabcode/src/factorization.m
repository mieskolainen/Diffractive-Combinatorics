%% Factorization of rapidity combinations
%
% mikael.mieskolainen@cern.ch, 2018
% Licensed under the MIT License <http://opensource.org/licenses/MIT>.

function [F,fact_labels] = factorization(X, N, labels)

F = [];

C = createCBM(N, true);

% Create label strings
fact_labels = {};

found = 1;
for i = 1:size(C,1)

    first_ind = i;

    % Find the pair
    c = C(i,:);
    flip_c = fliplr(c);

    % 1. Find the mirror pair
    pair_ind = -1;
    for k = 1:size(C,1)
        if (sum(C(k,:) == flip_c) == N)
            pair_ind = k;
            break;
        end
    end

    % 2. Find the factorization triple
    triple_ind = -1;
    for k = 1:size(C,1)
        if ( norm((c + flip_c) - C(k, :)) == 0)
           triple_ind = k;
           break;
        end
    end
    [first_ind pair_ind triple_ind]

    % Fix the indexing
    first_ind  = first_ind - 1;
    pair_ind   = pair_ind - 1;
    triple_ind = triple_ind - 1;

    % Found a valid combination
    if (triple_ind > 0)
        for j = 1:size(X,2) % Loop over thresholds
            F(found,j) =  X(triple_ind,j) / sqrt( X(first_ind,j) * X(pair_ind,j) );
        end
        fact_labels{end+1} = sprintf(' %s', labels{triple_ind});
        found = found + 1;
    end
end

end