% Diffractive Cross Section extraction plots
%
% mikael.mieskolainen@cern.ch, 2018
% Licensed under the MIT License <http://opensource.org/licenses/MIT>.

clear; close all;
addpath /home/user/cernbox/#matlabcodes
addpath ./src

% Model selection
mc = {'Pythia-6_(17h7a)', 'Phojet_(17h7b)'};
mc_model = 1;

% ========================================================================
% One chooses here manually the Pomeron intercept refit value for all
% values pre-calculated by C++ code.
POMERON_DELTA = 0.09;

% Numerical treshold to find out the right rows of .csv file. Check
% produced plots visually such that right indices were picked!
THRESHOLD = 0.025;
% ========================================================================

% Define indices of .csv files produced by C++ code
xs_ind = 2;              % Cross Section
xs_stat_ind = 3;         % Statistical uncertainty
xs_lumi_ind = 4;         % Luminosity uncertainty
eff_ind = 5;             % Fiducial acceptance
logL_ind = 6;            % -log(likelihood)
KLdiv_ind = 7;           % Kullback-Leibler
KSerr_ind = 8;           % Kolmogorov-Smirnov
chi2_ind  = 9;           % chi^2
deltaPomeron_ind = 10;   % Pomeron intercept Delta_P
deltaY_ind = 11;         % Minimum rapidity gap

% Runlist
runID = [274593 274594 274595];

for run = runID
for level = 1:3

fig1 = figure;

% Detector level
filename = sprintf('../figures_xsec/%d/CrossSections/Extraction_level_%d_Input_Data-%d_Model_%s.csv', run, level, run, mc{mc_model});
X = csvread(filename, 1,0);

CS = [];
CS_stat = [];
CS_lumi = [];

% Loop over processes
for id = [0 1 2 4] % sdl,sdr,dd,nd
    
    x = X(X(:,1) == id, :);
    
    % Find the right indices
    ind = (abs(x(:,deltaPomeron_ind) - POMERON_DELTA)/POMERON_DELTA < THRESHOLD);
    
    xval = x(ind, deltaY_ind);
    
    if (level == 1)
        eff = x(ind, eff_ind);
    end
    if (level == 2)
       eff = 1; 
    end
    cross_section      = x(ind, xs_ind) ./ eff ;
    cross_section_stat = x(ind, xs_stat_ind) ./ eff ;
    cross_section_lumi = x(ind, xs_lumi_ind) ./ eff ;
    
    CS      = [CS; cross_section(:)'];
    CS_stat = [CS_stat; cross_section_stat(:)'];
    CS_lumi = [CS_lumi; cross_section_lumi(:)'];
end

% Custom colors
colors = [1 0 0;
          0.2 0.8 0;
          0 0 0.8;
          0.4 0.4 0.4];

h = {};
transparency = 0.1;
for i = 1:size(CS,1)

    if (i == 4) % ND process
        f = 0.1; % Non-diffractive with / 10 for visualization
    else
        f = 1.0;
    end
    
    h{i} = plot(xval, CS(i,:)*f, 'color', colors(i,:), 'linewidth', 1); hold on;
    plotfill(xval, (CS_stat(i,:) + CS(i,:))*f, (CS(i,:) - CS_stat(i,:))*f, colors(i,:)*0.001, [1 1 1], transparency*1.1);
    plotfill(xval, (CS_lumi(i,:) + CS(i,:))*f, (CS(i,:) - CS_lumi(i,:))*f, colors(i,:), [1 1 1], transparency*0.5);
end

axis square;

index = 1; % We take the first (<DeltaY> = 0), no large difference
sigma_inel     = sum(CS(:,index));
sigma_inel_err = norm([sum(CS_stat(:,index)); sum(CS_lumi(:,index))]);

l = legend([h{1:4}], {'SD$_L$','SD$_R$','DD','ND / 10'});
set(l,'interpreter','latex'); legend('boxoff');
xlabel('$\langle \Delta Y \rangle_{\min}$ cutoff definition','interpreter','latex');
if (level == 1)
    ylabel('Total Cross Section (mb)','interpreter','latex');
    title(sprintf('$\\sqrt{s} = 13$ TeV, $\\, \\Delta_P^{fit} = %0.3f$, $\\, \\sigma_{inel}^{tot} = %0.1f \\pm %0.1f$ mb', POMERON_DELTA, sigma_inel, sigma_inel_err),'interpreter','latex');
end
if (level == 2)
    ylabel('Fiducial Cross Section (mb)','interpreter','latex');
    title(sprintf('$\\sqrt{s} = 13$ TeV, $\\, \\Delta_P^{fit} = %0.3f$, $\\, \\sigma_{inel}^{fid} = %0.1f \\pm %0.1f$ mb', POMERON_DELTA, sigma_inel, sigma_inel_err),'interpreter','latex');
end

set(gca,'XTick', 0:0.5:max(xval));

if (level == 1)
    axis([0 inf 4.5 9]);
end
if (level == 2)
    axis([0 inf 2.5 8]);
end

outputfile = sprintf('Run_%d_cross_sections_level_%d', run, level);
print(fig1, sprintf('./xsfigs/%s.pdf', outputfile), '-dpdf');
system(sprintf('pdfcrop --margins 10 ./xsfigs/%s.pdf ./xsfigs/%s.pdf', outputfile, outputfile));

end
end

