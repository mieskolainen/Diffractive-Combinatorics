% Toy test of Mass Distribution re-weighting
% mikael.mieskolainen@cern.ch, 2018
% Licensed under the MIT License <http://opensource.org/licenses/MIT>.

clear;
close all;

% Invariant Mass Squared
M2 = linspace(1.08^2, 10^2, 10000);

% Original distribution
f_pythia = @(M2, delta) 1./((M2).^(1 + delta));

% Reweighted distribution
f_new = @(M2, delta) 1./((M2).^(1 + delta*2));

% Plot orignal
delta = 0.08;
figure;
plot(sqrt(M2), f_pythia(M2, delta), 'b'); axis square;
set(gca,'yscale','log');
set(gca,'xscale','log');
xlabel('$M$ (GeV)', 'interpreter', 'latex');
ylabel('$d\\\sigma/dM$ (mb/GeV)', 'interpreter', 'latex');

hold on;
loglog(sqrt(M2), f_pythia(M2, delta*2), 'b--');

hold on;
loglog(sqrt(M2), f_pythia(M2, delta/2), 'b--');

hold on;
loglog(sqrt(M2), f_pythia(M2, delta) .* (M2).^(-delta), 'r--');

