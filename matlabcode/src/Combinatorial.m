% - Plot combinatorial partial cross sections (PLOT_ON = true)
% - Fit total inelastic                       (PLOT_ON = false)
%
% mikael.mieskolainen@cern.ch, 2018
% Licensed under the MIT License <http://opensource.org/licenses/MIT>.

function [mins, acceptances, METRICS] = Combinatorial(param, PLOT_ON, DATA_ON)

% ========================================================================
% Set parameters
run          = param.run;
unfold_model = param.model;

% Unfolding iteration range
hih = param.hih;
cen = param.cen;
low = param.low;

relative_lumi_uncertainty = param.relative_lumi_uncertainty;
sigma_inel_fid            = param.sigma_inel_fid;

basepath = param.basepath;
% ========================================================================

% This is beam-gas substracted and unfolded rates -> particle level   
for unfold_iter = 1:25
    str = sprintf('%s/%d/Ascii/Data-%d_x_unfolded_rates_iter_%d_model_%d.csv', ...
                            basepath, run, run, unfold_iter, unfold_model);
    ASCII = csvread(str, 1,0);
    DATA(:,unfold_iter) = ASCII(:,2);
end

%
DATA(1,:) = 0; % Null the extrapolation

% Dimension
N = log2(size(DATA,1));

% Combination legends
ylabels = makelegend(N, true);


% Read in MC from HepMC2 processing
basepath   = '/home/user/cernbox/CRMC/analyzer/output/';
filenames  = {'epos_lhc.dat', 'pythia6.dat', 'qgsjetII.dat','pythia8230_mbr.dat','pythia8230_default.dat', 'phojet.dat','sibyll.dat'};
labelnames = {'Epos-LHC', 'Pythia-6.x', 'QGSJet-II-04', 'Pythia-8.23$_{MBR}$', 'Pythia-8.23', 'Phojet-1.2', 'Sibyll-2.1'};

XX = cell(length(filenames),1);

acceptances = zeros(length(filenames), 1);
for i = 1:length(filenames)
    
    % Read MC
    XX{i} = dlmread([basepath filenames{i}], ' ');
    
    % Calculate fiducial acceptance
    acceptances(i) = sum(XX{i}(2:end,1)) / sum(XX{i}(:,1) );
    
    fprintf('MC %s : \tTotal events = %d \n', labelnames{i}, sum(XX{i}(:,1)) );
end

% Collect different unfolded results from C++ code
SIGMA     = zeros(size(DATA));
SIGMA_ERR = zeros(size(DATA));

for i = 1:25
    X = DATA(:,i); % Different unfolding iteration
    SIGMA(:,i)     = X        / sum(X) * sigma_inel_fid;
    SIGMA_ERR(:,i) = sqrt(X)  / sum(X) * sigma_inel_fid;
end


% Choose central point (based on results of AnalysisRegularization)

% Indices set here
high_ind    = 1;
central_ind = 2;
low_ind     = 3;

SIGMA_(:,high_ind)    = SIGMA(:,hih);
SIGMA_(:,central_ind) = SIGMA(:,cen);
SIGMA_(:,low_ind)     = SIGMA(:,low);

P_ = SIGMA_(2:end, central_ind);
P_ = P_ / sum(P_);

% Plot combinations in 4 different figures
ind_list  = [ 1 16;
             17 32;
             33 48;
             49 64];

if (PLOT_ON)
    MC_cross_section = param.MC_cross_section; % NOTE! TOTAL INELASTIC
else
    MC_cross_section = linspace(70, 74, 100);  % NOTE! FIDUCIAL INELASTIC
end

METRICS = zeros(length(XX), 5, length(MC_cross_section));

% Loop over MC total inelastic value
for ss = 1:length(MC_cross_section)

sigma_inel_fid_MC = MC_cross_section(ss);

SE = cell(7,1);
S2 = cell(7,1);
KL = cell(7,1);

% Loop over data in blocks of indices
for i = 1:size(ind_list,1)
    
    if (PLOT_ON)
    fig = figure('units','normalized','outerposition',[0 0 1 1]);
    end
    
    % Indices in reverse order
    ind =  ind_list(i,2):-1:ind_list(i,1);
    
    % Choose indices for this subfigure (1 out of 4)
    yval = fliplr(ind);
    
    % Draw unfolding + luminosity uncertainty, loop over all points
    systematic = zeros(length(yval),1);
    for k = 1:length(yval)
        
        color = [0.95 0.9 0.9];
        
        % Upper [central ... high)
        unfold_high = abs(SIGMA_(ind(k),high_ind) - SIGMA_(ind(k),central_ind));
        
        % Quadrature unfold (+) lumi
        systematic_high = norm([unfold_high SIGMA_(ind(k),central_ind)*relative_lumi_uncertainty]);
        
        if (DATA_ON)
        if (PLOT_ON)
        rectangle('Position', [SIGMA_(ind(k),central_ind) yval(k)-0.5 systematic_high 1.0], ...
            'FaceColor', color, 'EdgeColor', color, 'linewidth', 0.1);
        end
        end
        % Lower (lower ... central]
        unfold_low = abs(SIGMA_(ind(k),low_ind) - SIGMA_(ind(k),central_ind));
        
        % Quadrature unfold (+) lumi
        systematic_low = norm([unfold_low SIGMA_(ind(k),central_ind)*relative_lumi_uncertainty]);
        
        if (DATA_ON)
        if (PLOT_ON)
        rectangle('Position', [SIGMA_(ind(k),central_ind)-systematic_low yval(k)-0.5 systematic_low 1.0], ...
            'FaceColor', color, 'EdgeColor', color, 'linewidth', 0.1);
        end
        end
        
        % Take maximum as the total systematic
        systematic(k) = max([systematic_high, systematic_low]);
        
    % Draw statistical uncertainty, loop over all points
    if (DATA_ON)
    if (PLOT_ON)
        hold on;    
        for k = 1:length(yval)

            color = ones(3,1)*0.8;
            width = SIGMA_ERR(ind(k),cen);

            % Upper [central ... high)
            rectangle('Position', [SIGMA_(ind(k),central_ind) yval(k)-0.5 width 1.0], ...
                'FaceColor', color, 'EdgeColor', color, 'linewidth', 0.1);

            % Lower (lower ... central]
            rectangle('Position', [SIGMA_(ind(k),central_ind)-width yval(k)-0.5 width 1.0], ...
                'FaceColor', color, 'EdgeColor', color, 'linewidth', 0.1);
        end
    end
    end
       
    end
    %% Monte Carlo values
    
    % Loop over MC generators
    h = {};
    COLORMAT = get(gca, 'ColorOrder');
    for f = 1:length(XX)
        
        % Take first column
        MC = XX{f}(:,1);
        
        % Null the first bin (UNVISIBLE, NOT USED IN THE FIT!)
        MC(1) = 0;
        
        % Binomial (multinomial) counting error
        n = sum(MC);
        phat = MC ./ n;
        MC_err = sqrt( phat.*(1-phat) / n);
        
        % Scale to cross sections
        MC_ = MC / n * sigma_inel_fid_MC;
        MC_err_ = MC_err * sigma_inel_fid_MC;
        
        % Residual
        residual =  SIGMA_(ind, central_ind) - MC_(ind);
        SE{f} = [SE{f}; flipud(residual)];
        
        % Residual/Uncertainty
        alpha = (sigma_inel_fid - sigma_inel_fid_MC)/sigma_inel_fid;
        delta = relative_lumi_uncertainty;
        
        % Chi2 residuals
        chi2res = ( (SIGMA_(ind, central_ind) - MC_(ind))./ sqrt(SIGMA_ERR(ind,cen).^2 + systematic.^2) ).^2 + ...
                + alpha^2/delta^2;
        S2{f} = [S2{f}; flipud(chi2res)];
        
        % Kullback-Leibler
        MCP = MC(2:end); MCP = MCP / sum(MCP);
        
        if (i == 1)
            klval = -P_ .* log2(MCP ./ P_);
            KL{f} = [KL{f}; flipud(klval)];
        end
        
        % Plot
        if (PLOT_ON)
        h{f} = plot(MC_(ind), yval, 's', 'Markersize', 9, 'Color', COLORMAT(f,:)); hold on;
        errorbar(MC_(ind), yval, zeros(size(yval)), zeros(size(yval)), MC_err_(ind), MC_err_(ind), ...
            's', 'Markersize', 0.1, 'CapSize', 4, 'Color', COLORMAT(f,:));
        end
    end
    
    %% Central unfolded value
    
    if (DATA_ON)
    if (PLOT_ON)
    for k = length(yval):-1:1
        
        color = [0 0 0];
        
        % Unvisible estimate with red (combination 0)
        if (i == 1 && k == 16) 
            color = [1 0 0];
        end
        hc = errorbar(SIGMA_(ind(k), central_ind), yval(k), 0.5, 0.5, ...
            's', 'Color', color, 'MarkerFaceColor', color, 'CapSize', 0); hold on;
    end
    end
    end
    
    %% Axis Ticks
    
    if (PLOT_ON)
    % Y-axis
    yticks(fliplr(ind));
    yticklabels(ylabels(ind));
    
    % X-axis
    xlabel('$\sigma_k$ (mb)','interpreter','latex');
    xticks([1e-3 1e-2 1e-1 1 1e1 1e2]);
    ax = gca;
    ax.TickLabelInterpreter = 'latex';
    
    % Plot grid
    grid on;
    %grid minor;
    
    % Grid on
    set(gca,'xscale','log');
    axis([1e-3 100 ind(end)-1 ind(1)+1]);
    
    % Plot aspect ratio
    pbaspect([1 2.5 1]);
    
    % Legend
    if (DATA_ON)
        if (i == 1 || i == 3)
            l = legend([hc h{1} h{2} h{3} h{4} h{5} h{6} h{7}], ...
                'ALICE 13 TeV', labelnames{1}, labelnames{2}, labelnames{3}, labelnames{4}, labelnames{5}, labelnames{6}, labelnames{7} );
            set(l,'interpreter','latex');
            set(l,'Position',[0.513 0.209 0.082 0.138]);
            %set(l,'Position',[0.513 0.688 0.082 0.138]); 
        end
    else
        if (i == 1 || i == 3)
        l = legend([h{1} h{2} h{3} h{4} h{5} h{6} h{7}], ...
            labelnames{1}, labelnames{2}, labelnames{3}, labelnames{4}, labelnames{5}, labelnames{6}, labelnames{7} );
        set(l,'interpreter','latex');
        set(l,'Position',[0.513 0.209 0.082 0.138]);
        %set(l,'Position',[0.513 0.688 0.082 0.138]); 
        end
    end
    
    % Print
    outputfile = sprintf('comb%d.pdf', i);
    print(fig, sprintf('./combfigs/%s', outputfile), '-dpdf');
    system(sprintf('pdfcrop --margins 10 ./combfigs/%s ./combfigs/%s', outputfile, outputfile));
    end
end

% Merge pdfs
if (PLOT_ON)
system('gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=./combfigs/COMBINED.pdf ./combfigs/comb1.pdf ./combfigs/comb2.pdf ./combfigs/comb3.pdf ./combfigs/comb4.pdf');
end

%% Calculate metrics
CHI2 = zeros(length(SE),1);
CABS = CHI2;
CTUK = CHI2;
CHUB = CHI2;
CKL  = CHI2;

% The set of combinations used to calculate metric
indexset = 2:2^N;

metriclabels = {'$\chi^2$', 'MAD', 'Tukey bi-weight', 'Huber weight', 'KL'};
format long;

% Calculate cost values
for f = 1:length(SE)
    
    % Chi^2
    CHI2(f) = sum(S2{f}(indexset));   
    
    % Absolute
    CABS(f) = madcost(SE{f}(indexset));
    
    % Tukey
    CTUK(f) = tukeycost(S2{f}(indexset));
    
    % Huber
    CHUB(f) = hubercost(S2{f}(indexset));
    
    % Kullback-Leibler
    KL{f}(isinf(KL{f})) = 0;
    CKL(f)    = sum(KL{f});
end

fprintf('MC(sigma_inel_fid) = %0.2f mb \n\n', sigma_inel_fid_MC);
fprintf('chi2 abs tukey huber KL \n');
disp([CHI2 CABS CTUK CHUB CKL]);
METRICS(:,:,ss) = [CHI2 CABS CTUK CHUB CKL];

end

mins = [];

% Plot the behavior
if (~PLOT_ON)
for i = 1:size(METRICS,2) % Loop over metrics
    
    fig1 = figure('units','normalized','outerposition',[0 0 0.6 0.6]);
    mins = zeros(size(METRICS,1), 1);
    
    % Loop over MC
    for k = 1:size(METRICS,1)
        subplot(2,4,k);
        
        % Plot
        costval = squeeze(METRICS(k,i,:));
        xvalues = MC_cross_section / acceptances(k);
        plot(xvalues, costval);
        
        [a,b] = min(costval);
        mins(k) = MC_cross_section(b);
        axis square;
        title(sprintf('%s', labelnames{k}),'interpreter','latex');
        
        ylabel(sprintf('%s', metriclabels{i}),'interpreter','latex');
        xlabel('$\hat{\sigma}_{inel}^{extrapolated}$ (mb)','interpreter','latex');
        set(gca,'yscale','log');
        axis([min(xvalues) max(xvalues) -inf inf]);
        xticks(round(min(xvalues):1:max(xvalues),0));
        axis tight;
        hold on;
    end
    
    outputfile = sprintf('run_%d_total_fit_%d.pdf', run, i);
    print(fig1, sprintf('./combfigs/%s', outputfile), '-dpdf');
    system(sprintf('pdfcrop --margins 10 ./combfigs/%s ./combfigs/%s', outputfile, outputfile));
    
end
end
end



