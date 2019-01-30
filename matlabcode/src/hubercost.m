% Huber robust cost
%
% Input:  residual = difference vector (Nx1)
%
% mikael.mieskolainen@cern.ch, 2018
% Licensed under the MIT License <http://opensource.org/licenses/MIT>.

function cost = hubercost(residual)

% Asymptotic 95% CL efficiency of the estimate (see Robust Regression theory)
c = 1.345;

cost = 0;
for i = 1:length(residual)

    if (abs(residual(i)) <= c)
        cost = cost + 0.5*residual(i)^2;
    else
        cost = cost + c * (abs(residual(i) - c/2));
    end
end

end