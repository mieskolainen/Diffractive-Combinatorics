% Diffractive Cross Section extraction plots
%
% mikael.mieskolainen@cern.ch, 2019
% Licensed under the MIT License <http://opensource.org/licenses/MIT>.

clear; close all;
addpath /home/user/cernbox/#matlabcodes
addpath ./src

% Model selection
mc = {'Pythia-6_(17h7a)', 'Phojet_(17h7b)'};

for mc_model = 1:length(mc)

% ========================================================================
% One chooses here manually the Pomeron intercept refit value for all
% values pre-calculated by C++ code.
POMERON_DELTA = 0.09;

% Numerical treshold to find out the right rows of .csv file. Check
% produced plots visually such that right indices were picked!
THRESHOLD = 0.025;
% ========================================================================

% Define indices of .csv files produced by C++ code
xs_ind      = 2;         % Cross Section
xs_stat_ind = 3;         % Statistical uncertainty
xs_lumi_ind = 4;         % Luminosity uncertainty
eff_ind     = 5;         % Fiducial acceptance
logL_ind    = 6;         % -log(likelihood)
KLdiv_ind   = 7;         % Kullback-Leibler
KSerr_ind   = 8;         % Kolmogorov-Smirnov
chi2_ind    = 9;         % chi^2
deltaPomeron_ind = 10;   % Pomeron intercept Delta_P
deltaY_ind = 11;         % Minimum rapidity gap

% Runlist [the first defines the central point, others give run-by-run variations]
runID = [274595 274593 274594];

% Extraction level
for level = 1:2

CS_all      = {};
CS_stat_all = {};
CS_lumi_all = {};

for run = runID

    % Detector level
    filename = sprintf('../figures_xsec/%d/CrossSections/Extraction_level_%d_Input_Data-%d_Model_%s.csv', ...
                       run, level, run, mc{mc_model});
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
        
        fprintf('xs[%d] : Efficiency = \n', id);
        disp(eff);
        
        cross_section      = x(ind, xs_ind) ./ eff ;
        cross_section_stat = x(ind, xs_stat_ind) ./ eff ;
        cross_section_lumi = x(ind, xs_lumi_ind) ./ eff ;
        
        CS      = [CS; cross_section(:)'];
        CS_stat = [CS_stat; cross_section_stat(:)'];
        CS_lumi = [CS_lumi; cross_section_lumi(:)'];
    end
    
    CS_all{end+1}      = CS;
    CS_stat_all{end+1} = CS_stat;
    CS_lumi_all{end+1} = CS_lumi;
end % Over run

%%
% Custom colors
colors = [1 0 0;
          0.2 0.8 0;
          0 0 0.8;
          0.4 0.4 0.4];

h = {};
transparency = 0.1;

% Choose the first run, others give systematic variations
CS      = CS_all{1};
CS_stat = CS_stat_all{1};
CS_lumi = CS_lumi_all{1};

% Run-by-Run variations
CS_rbr_min = zeros(size(CS_all{1}));
CS_rbr_max = zeros(size(CS_all{1}));

for ID1 = 1
    for ID2 = 1:length(runID)
        [CS_rbr_min, CS_rbr_max] = minmaxvar(CS_all{ID1}, CS_all{ID2}, CS_rbr_min, CS_rbr_max);
    end
end


%% Inelastic

% We take the first (<DeltaY> = 0), no large difference
index = 1;
sigma_inel     = sum(CS_all{1}(:,index));

% Run-by-Run
sigma_inel_rbr = sum( (abs(CS_rbr_min(:,index)) + abs(CS_rbr_min(:,index))) / 2 ); % avg between max/min

% Error in quadrature
sigma_inel_err = norm([sum(CS_stat_all{1}(:,index)); sum(CS_lumi_all{1}(:,index)); sigma_inel_rbr]);


%% Plotting

fig1 = figure;

for i = 1:size(CS,1)
    
    if (i == 4)  % ND process
        f = 0.1; % Non-diffractive with / 10 for visualization
    else
        f = 1.0;
    end
    
    h{i} = plot(xval, CS(i,:)*f, 'color', colors(i,:), 'linewidth', 1); hold on;
    plotfill(xval, (CS(i,:) + CS_stat(i,:))*f,    (CS(i,:) - CS_stat(i,:))*f,    colors(i,:)*0.001, [1 1 1], transparency*1.1);
    plotfill(xval, (CS(i,:) + CS_lumi(i,:))*f,    (CS(i,:) - CS_lumi(i,:))*f,    colors(i,:), [1 1 1],       transparency*0.5);
    plotfill(xval, (CS(i,:) + CS_rbr_max(i,:))*f, (CS(i,:) + CS_rbr_min(i,:))*f, colors(i,:), [1 1 1],       transparency*0.5);
end

axis square;

l = legend([h{1:4}], {'SD$_L$','SD$_R$','DD','ND / 10'});
set(l,'interpreter','latex'); legend('boxoff');
xlabel('$\langle \Delta Y \rangle_{\min}$ cutoff definition','interpreter','latex');
if (level == 1)
    ylabel('Total Cross Section (mb)','interpreter','latex');
    title(sprintf('$\\sqrt{s} = 13$ TeV, $\\, \\Delta_P^{fit} = %0.3f$, $\\, \\sigma_{inel}^{tot} = %0.1f \\pm %0.1f$ mb', ...
        POMERON_DELTA, sigma_inel, sigma_inel_err),'interpreter','latex');
end
if (level == 2)
    ylabel('Fiducial Cross Section (mb)','interpreter','latex');
    title(sprintf('$\\sqrt{s} = 13$ TeV, $\\, \\Delta_P^{fit} = %0.3f$, $\\, \\sigma_{inel}^{fid} = %0.1f \\pm %0.1f$ mb', ...
        POMERON_DELTA, sigma_inel, sigma_inel_err),'interpreter','latex');
end

set(gca,'XTick', 0:0.5:max(xval));

if (level == 1)
    axis([0 inf 3 9]);
end
if (level == 2)
    axis([0 inf 3 9]);
end
set(gca,'YTick', 3:0.5:9);

outputfile = sprintf('Run_%d_cross_sections_level_%d_mcmodel_%d', runID(1), level, mc_model);
print(fig1, sprintf('./xsfigs/%s.pdf', outputfile), '-dpdf');
system(sprintf('pdfcrop --margins 10 ./xsfigs/%s.pdf ./xsfigs/%s.pdf', outputfile, outputfile));

end
end

%%
system('source copyplots.sh');


