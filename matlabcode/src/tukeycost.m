% Tukey biweight robust cost
%
% Input:  residual = difference vector (Nx1)
%
% mikael.mieskolainen@cern.ch, 2018
% Licensed under the MIT License <http://opensource.org/licenses/MIT>.

function cost = tukeycost(residual)

% Asymptotic 95% CL efficiency of the estimate (see Robust Regression theory)
c = 4.685; 

cost = 0;
for i = 1:length(residual)

    if (abs(residual(i)) <= c)
        cost = cost + ( 1 - (1 - (residual(i)/c)^2)^3 );
    else
        cost = cost + 1;
    end
end

end