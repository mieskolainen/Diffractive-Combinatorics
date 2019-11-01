% Plot Partial Combinatorial Cross Sections
%
% mikael.mieskolainen@cern.ch, 2019
% Licensed under the MIT License <http://opensource.org/licenses/MIT>.

clear; close all;
addpath ./src
param.basepath = '../figures_xsec';

% Unfolding iteration range
param.hih = 6;           % Maximum iteration
param.cen = 4;           % Central value
param.low = 2;           % Minimum iteration

% *** SET HERE MANUALLY BASED ON C++ CODE + EXTERNAL INFO ***
param.relative_lumi_uncertainty = 1.2 / 57.8; % V0-AND vdM +- / V0-AND vdM value
param.sigma_inel_fid = 72;                    % Measured total fiducial in data
param.sigma_inel_MC  = 80;                    % Total inelastic in MC

% First one is the central, others are used for run-by-run variations
param.runs = [274595, 274593, 274594];

for mode = [true false]

PLOT_ON = mode; % if false, only table printed out in proper ID order
DATA_ON = true;
FIT_ON  = true;
for model = 1
    param.model = model;
    [mins, acceptances, METRICS] = Combinatorial(param, PLOT_ON, DATA_ON, FIT_ON);
end


%% Print out inelastic box figure

if (FIT_ON)

fig = figure;

mu_fid = param.sigma_inel_fid;
uncert_fid = param.relative_lumi_uncertainty * param.sigma_inel_fid;

color = [0.95 0.9 0.9];

% Upper [central ... high)
rectangle('Position', [mu_fid 0.67 uncert_fid 0.16], ...
'FaceColor', color, 'EdgeColor', color, 'linewidth', 0.1);

% Lower (lower ... central]
rectangle('Position', [mu_fid-uncert_fid 0.67 uncert_fid 0.16], ...
'FaceColor', color, 'EdgeColor', color, 'linewidth', 0.1); hold on;

errorbar(mu_fid, 0.75, 0, 0, 0.5, 0.5,'ks', 'MarkerSize', 11); hold on;

text(mu_fid*0.993,  0.88, '$\sigma_{\mathrm{fid}}^{measured}$','interpreter','latex','fontsize',22);
text(mu_fid*0.978,  0.60, 'Fiducial (approx):','interpreter','latex','fontsize', 12);
text(mu_fid*0.978,  0.55, '$-7.0 < \eta_{ch} < 6.3$','interpreter','latex','fontsize', 12);
text(mu_fid*0.978,  0.51, '$p_T > 50$ MeV','interpreter','latex','fontsize', 12);

text(mu_fid*0.987,  0.69, 'Luminosity','interpreter','latex','fontsize',10);
plot(mu_fid,        0.75, 'ks', 'MarkerSize', 12); hold on;

color = ones(3,1)*0.9;

% Robust estimates
MC_estimate = mins ./ acceptances;

mu_tot     = 78.4; % From Pythia-6 re-tune
tot_uncert = median(abs(MC_estimate - median(MC_estimate))) * 1.4826;

% Upper [central ... high)
%rectangle('Position', [mu_tot 0.67 tot_uncert 0.16], 'FaceColor', color, 'EdgeColor', color, 'linewidth', 0.1);

% Lower (lower ... central]
%rectangle('Position', [mu_tot-tot_uncert 0.67 tot_uncert 0.16], 'FaceColor', color, 'EdgeColor', color, 'linewidth', 0.1);

errorbar(mu_tot, 0.75, 0, 0, mu_tot-min(MC_estimate), max(MC_estimate)-mu_tot, 'ko', 'MarkerSize', 11);

text(mu_tot*0.99, 0.88,  '$\hat{\sigma}_{\mathrm{inel}}^{extrapolation}$','interpreter','latex','fontsize',22);
text(mu_tot*0.97, 0.7,   'Max Likelihood Pythia 6 Re-Tune','interpreter','latex','fontsize', 9.4);

text(min(MC_estimate), 0.78,  'MC-min','interpreter','latex','fontsize', 8);
text(max(MC_estimate), 0.78, ' MC-max','interpreter','latex','fontsize', 8);

text(77.5, 1.035, '$\sqrt{s} = 13$ TeV', 'interpreter', 'latex', 'fontsize', 18);

xlabel('$\sigma$ (mb)', 'interpreter', 'latex');

axis([69 82 0.4 1.1]);
xticks([69:82]);

% Remove y-axis
set(gca,'YTickLabel',{' '})

% Border
box on;

outputfile = 'inelastic.pdf';
print(fig, sprintf('./combfigs/%s', outputfile), '-dpdf');
system(sprintf('pdfcrop --margins 10 ./combfigs/%s ./combfigs/%s', outputfile, outputfile));

end

end
%%
system('source copyplots.sh');
