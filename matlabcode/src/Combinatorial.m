% ------------------------------------------------------------------------
%
% Plot combinatorial partial cross sections (PLOT_ON = true)
% Fit total inelastic                       (FIT_ON  = true and DATA_ON = true)
% 
% mikael.mieskolainen@cern.ch, 2019
% Licensed under the MIT License <http://opensource.org/licenses/MIT>.

function [minsxs, acceptances, METRICS] = Combinatorial(param, PLOT_ON, DATA_ON, FIT_ON)

fprintf('\n');

minsxs = [];
acceptances = [];
METRICS = [];

if (FIT_ON && ~DATA_ON)
    fprintf('Combinatorial: Error with FIT_ON = true but DATA_ON = false \n');
    return;
end

% Read in MC from HepMC2 processing
mcbasepath   = '/home/user/cernbox/CRMC/analyzer/output/';
filenames  = {'epos_lhc.dat', 'pythia6.dat', 'qgsjetII.dat','pythia8230_mbr.dat','pythia8230_default.dat', 'phojet.dat','sibyll.dat'};
labelnames = {'Epos-LHC', 'Pythia-6.x', 'QGSJet-II-04', 'Pythia-8.23$_{MBR}$', 'Pythia-8.23', 'Phojet-1.2', 'Sibyll-2.1'};

% How many event generators
NUMBER_OF_MC = length(labelnames);
EVENTS_MC    = cell(length(filenames),1);
acceptances = zeros(length(filenames), 1);

for i = 1:length(filenames)
    
    % Read MC
    EVENTS_MC{i} = dlmread([mcbasepath filenames{i}], ' ');
    
    % Calculate fiducial acceptance
    acceptances(i) = sum(EVENTS_MC{i}(2:end,1)) / sum(EVENTS_MC{i}(:,1) );
    
    fprintf('MC %s : \tTotal events = %d [Fiducial acceptance = %0.2f] \n', labelnames{i}, sum(EVENTS_MC{i}(:,1)), acceptances(i) );
end

% Dimension
N = log2(size(EVENTS_MC{1}, 1));

% Combination legends
ylabels = makelegend(N, true);

% ------------------------------------------------------------------------
if (FIT_ON)
    % % NOTE! FIDUCIAL INELASTIC, put here visually suitable range
    MC_fid_inel_xsection = linspace(70, 74, 100);
    
    SE = cell(NUMBER_OF_MC, length(MC_fid_inel_xsection));
    S2 = cell(NUMBER_OF_MC, length(MC_fid_inel_xsection));
end

% ========================================================================
if (DATA_ON)

% Set parameters
unfold_model = param.model;

relative_lumi_uncertainty = param.relative_lumi_uncertainty;
sigma_inel_fid            = param.sigma_inel_fid;

basepath = param.basepath;
% ========================================================================
    
% Number of regularization iterations
REGITER = 25;

% This is beam-gas subtracted and unfolded rates -> particle level
EVENTS = {};

for ID = 1:length(param.runs)
    for unfold_iter = 1:REGITER
        str = sprintf('%s/%d/Ascii/Data-%d_x_unfolded_rates_iter_%d_model_%d.csv', ...
                                basepath, param.runs(ID), param.runs(ID), unfold_iter, unfold_model);
        ASCII = csvread(str, 1,0);
        EVENTS{ID}(:,unfold_iter) = ASCII(:,2);
    end
end

%
for ID = 1:length(EVENTS)
    EVENTS{ID}(1,:) = 0; % Null the extrapolation
end    

% Collect different unfolded results from C++ code
SIGMA     = zeros(size(EVENTS{1}));
SIGMA_ERR = zeros(size(EVENTS{1}));

% Over different unfolding iterations
for i = 1:REGITER
    X = EVENTS{1}(:,i);
    SIGMA(:,i)     = X       / sum(X) * sigma_inel_fid;
    SIGMA_ERR(:,i) = sqrt(X) / sum(X) * sigma_inel_fid;
end

% ------------------------------------------------------------------------
% Print out run-by-run
fprintf('RUN-by-RUN central values (events): \n');
for c = 1:2^N
   fprintf('%2d : ', c);
   for ID = 1:length(param.runs)
       fprintf('%0.2f \t', EVENTS{ID}(c, param.cen));
   end
   fprintf('\n');
end

% ------------------------------------------------------------------------
% Find out run-by-run relative minimum and maximum
RBR_relative_low  = ones(size(SIGMA));
RBR_relative_high = ones(size(SIGMA));

for ID1 = 1:length(param.runs)
    for ID2 = 1:length(param.runs)
        for c = 1:2^N
            for unfold_iter = 1:REGITER
                
                ratio = (EVENTS{ID2}(c, unfold_iter) / sum(EVENTS{ID2}(:, unfold_iter)) ) ...
                      / (EVENTS{ID1}(c, unfold_iter) / sum(EVENTS{ID1}(:, unfold_iter)) );
                
                if ratio > RBR_relative_high(c, unfold_iter)
                    RBR_relative_high(c, unfold_iter) = ratio;
                end
                if ratio < RBR_relative_low(c, unfold_iter)
                    RBR_relative_low(c, unfold_iter)  = ratio;
                end
            end
        end
    end
end

% Take into account the min/max uniform variations in terms of 1 sigma
% sqrt(12) comes from Var(uniform rnd) = (b-a)^2/12
RBR_relative_high = 1 + (RBR_relative_high - 1) / sqrt(12);
RBR_relative_low  = 1 - (1- RBR_relative_low)   / sqrt(12);

end % DATA_ON

% ========================================================================

% Plot combinations in 4 different figures
ind_list  = [ 1 16;
             17 32;
             33 48;
             49 64];

if (DATA_ON)
    fprintf('\n');
    fprintf('Run %d fiducial partial cross sections: \n', param.runs(1));
    
    % Open text output
    outputfile = sprintf('./combfigs/xstable_%d.tex', param.runs(1));
    fp = fopen(outputfile, 'w');
    
    fprintf('Writing .tex output to %s \n', outputfile);
end

% Loop over data in blocks of indices
for i = 1:size(ind_list,1)
    
    fprintf(fp,'\n');
    fprintf(fp,'\\begin{table}\n');
    fprintf(fp,'\\begin{center}\n');
    fprintf(fp,'\\renewcommand{\\arraystretch}{1.4}\n');
    fprintf(fp,'\\begin{tabular}{|cc|ccc|ccc|}\n');
    fprintf(fp,'\\hline\n');
    fprintf(fp,'& \\small{\\textsc{X-Section}} & value (mb) & stat & tot.syst & lumi & unfold & run-by-run \\\\ \n');
    fprintf(fp,'\\hline \n');

    if (PLOT_ON)
        fig = figure('units','normalized','outerposition',[0 0 1 1]);
    end
    
    if (PLOT_ON)
        % Indices in reverse order
        ind = ind_list(i,2):-1:ind_list(i,1);
    else
        % Normal order
        ind = ind_list(i,1):ind_list(i,2);
    end
    
    % Choose indices for this subfigure (1 out of 4)
    yval = fliplr(ind);
    
    if (DATA_ON)

        % Draw unfolding + luminosity uncertainty, loop over all points
        systematic = zeros(length(yval),1);
        for k = 1:length(yval)

            color = [0.95 0.9 0.9];
            
            % ============================================================
            % Collect uncertainties
            
            % Statistical uncertainty
            stat_err = SIGMA_ERR(ind(k) ,param.cen);

            % Luminosity uncertainty
            lumi_err = SIGMA(ind(k),     param.cen)*relative_lumi_uncertainty;
            
            % Unfolding uncertainty
            unfold_high = abs(SIGMA(ind(k), param.hih) - SIGMA(ind(k), param.cen));
            unfold_low  = abs(SIGMA(ind(k), param.low) - SIGMA(ind(k), param.cen));
            
            % Run-by-Run variations
            rbr_high    = abs( SIGMA(ind(k), param.cen)*RBR_relative_high(ind(k), param.cen) - SIGMA(ind(k), param.cen));
            rbr_low     = abs( SIGMA(ind(k), param.cen)*RBR_relative_low(ind(k),  param.cen) - SIGMA(ind(k), param.cen));
            
            % Total systematics in quadrature
            total_syst_high = sqrt(sum([lumi_err, unfold_high, rbr_high].^2));
            total_syst_low  = sqrt(sum([lumi_err,  unfold_low,  rbr_low].^2));
            
            % Take maximum as the total systematic
            systematic(k) = max([total_syst_low, total_syst_high]);
            
            % ============================================================
            
            % Print out
            fprintf(fp,'%d & %s & $%0.3f$ & $\\pm %0.3f$ & $_{-%0.3f}^{+%0.3f}$ & $\\pm %0.3f$ & $_{-%0.3f}^{+%0.3f}$ & $_{-%0.3f}^{+%0.3f}$ \\\\ \n', ...
                ind(k) - 1, ...
                ylabels{ind(k)}, ...
                SIGMA(ind(k), param.cen), ...
                stat_err, ...
                total_syst_low, total_syst_high, ...
                lumi_err, ...
                unfold_low, unfold_high, ...
                rbr_low, rbr_high);
            
            if (PLOT_ON)
                
                % Draw systematic uncertainty
                rectangle('Position', [SIGMA(ind(k), param.cen) yval(k)-0.5 total_syst_high 1.0], ...
                    'FaceColor', color, 'EdgeColor', color, 'linewidth', 0.1); hold on;

                rectangle('Position', [SIGMA(ind(k), param.cen)-total_syst_low yval(k)-0.5 total_syst_low 1.0], ...
                    'FaceColor', color, 'EdgeColor', color, 'linewidth', 0.1); hold on;
        
                % Draw statistical uncertainty, loop over all points
                color = ones(3,1)*0.8;
                width = SIGMA_ERR(ind(k), param.cen);

                % Upper [central ... high)
                rectangle('Position', [SIGMA(ind(k), param.cen) yval(k)-0.5 width 1.0], ...
                    'FaceColor', color, 'EdgeColor', color, 'linewidth', 0.1);

                % Lower (lower ... central]
                rectangle('Position', [SIGMA(ind(k), param.cen)-width yval(k)-0.5 width 1.0], ...
                    'FaceColor', color, 'EdgeColor', color, 'linewidth', 0.1);
            end
            
        end % Over combinations
    end % if DATA_ON
    
    
    %% Monte Carlo plotting
    
    if (PLOT_ON)
        
        % Plot first
        h = {};
        COLORMAT = get(gca, 'ColorOrder');

        % Loop over MC generators
        for f = 1:length(EVENTS_MC)

            % Take first column, correspond to minimum pt-threshold
            MC = EVENTS_MC{f}(:,1);
            
            % Binomial (multinomial) counting error
            n = sum(MC);
            phat = MC ./ n;
            MC_err = sqrt( phat.*(1-phat) / n);

            % Scale to cross sections with TOTAL inelastic values
            MC_ = MC / n * param.sigma_inel_MC;
            MC_err_ = MC_err * param.sigma_inel_MC;
            
            h{f} = plot(MC_(ind), yval, 's', 'Markersize', 9, 'Color', COLORMAT(f,:)); hold on;
            errorbar(MC_(ind), yval, zeros(size(yval)), zeros(size(yval)), MC_err_(ind), MC_err_(ind), ...
            's', 'Markersize', 0.1, 'CapSize', 4, 'Color', COLORMAT(f,:));

        end
    end
    
     %% Central unfolded value as the last marker
    
    if (DATA_ON && PLOT_ON)
        for k = length(yval):-1:1

            color = [0 0 0];

            % Unvisible estimate with red (combination 0)
            if (i == 1 && k == 16) 
                color = [1 0 0];
            end
            hc = errorbar(SIGMA(ind(k), param.cen), yval(k), 0.5, 0.5, ...
                's', 'Color', color, 'MarkerFaceColor', color, 'CapSize', 0); hold on;
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
   
    end % PLOT_ON

    
    %% FIT LOOP
    if (FIT_ON)
    
    % Loop over MC total inelastic value
    for ss = 1:length(MC_fid_inel_xsection)

        sigma_inel_fid_MC = MC_fid_inel_xsection(ss);
        
        %% Monte Carlo fit loop

        % Loop over MC generators
        for f = 1:length(EVENTS_MC)

            % Take first column
            MC       = EVENTS_MC{f}(:,1);

            % Null the first bin (UNVISIBLE, NOT USED IN THE FIT!)
            MC(1)    = 0;

            % Binomial (multinomial) counting error
            n        = sum(MC);
            %phat     = MC ./ n;
            %MC_err   = sqrt( phat.*(1-phat) / n);
            
            % Scale to cross sections with total fiducial values
            MC_      = MC / n * sigma_inel_fid_MC;
            %MC_err_  = MC_err * sigma_inel_fid_MC;

            % Residual
            residual =  SIGMA(ind, param.cen) - MC_(ind);
            
            if (PLOT_ON)
                residual = flipud(residual);
            end
            SE{f,ss}    = [SE{f,ss}; residual];

            % Residual/Uncertainty
            alpha    = (sigma_inel_fid - sigma_inel_fid_MC)/sigma_inel_fid;
            delta    = relative_lumi_uncertainty;

            % Chi2 residuals
            chi2res  = ( (SIGMA(ind, param.cen) - MC_(ind))./ sqrt(SIGMA_ERR(ind, param.cen).^2 + systematic.^2) ).^2 + ...
                    + alpha^2/delta^2;
            
            if (PLOT_ON)
                chi2res = flipud(chi2res);
            end
            S2{f,ss} = [S2{f,ss}; chi2res];
        end
        
    end % Loop over cross section values
    
    end % If FIT_ON
    
    fprintf(fp, '\\hline \n');
    fprintf(fp, '\\end{tabular} \n');
    fprintf(fp, '\\caption{Fiducial partial cross sections [%d,%d].} \n', ind_list(i,1)-1, ind_list(i,2)-1);
    fprintf(fp, '\\label{table:xstable%d} \n', i);
    fprintf(fp, '\\end{center} \n');
    fprintf(fp, '\\end{table} \n');
    
end % Over different combinations

if (PLOT_ON)
   system('gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=./combfigs/COMBINED.pdf ./combfigs/comb1.pdf ./combfigs/comb2.pdf ./combfigs/comb3.pdf ./combfigs/comb4.pdf');
end


%% FITS
if (FIT_ON)
    
NUMBER_OF_METRICS = 4;
METRICS = zeros(NUMBER_OF_MC, length(MC_fid_inel_xsection), NUMBER_OF_METRICS);

% The set of combinations used to calculate metric
indexset = 2:2^N;

metriclabels = {'$\chi^2$', 'MAD', 'Tukey bi-weight', 'Huber weight'};
format long;

% Calculate metrics
CHI2 = zeros(NUMBER_OF_MC, length(MC_fid_inel_xsection));
CABS = CHI2;
CTUK = CHI2;
CHUB = CHI2;

% Over each fiducial inelastic hypothesis
for ss = 1:length(MC_fid_inel_xsection)

    % Get costs for each MC
    for f = 1:size(SE,1)

        % Chi^2
        CHI2(f,ss) =       sum(S2{f,ss}(indexset));   

        % Absolute
        CABS(f,ss) =   madcost(SE{f,ss}(indexset));

        % Tukey
        CTUK(f,ss) = tukeycost(S2{f,ss}(indexset));

        % Huber
        CHUB(f,ss) = hubercost(S2{f,ss}(indexset));
    end
end

METRICS(:,:,1) = CHI2;
METRICS(:,:,2) = CABS;
METRICS(:,:,3) = CTUK;
METRICS(:,:,4) = CHUB;

% Plot the behavior
for i = 1:NUMBER_OF_METRICS
    
    fig1 = figure('units','normalized','outerposition',[0 0 0.6 0.6]);
    minsxs = zeros(size(METRICS,1), 1);
    
    % Loop over MC
    for k = 1:NUMBER_OF_MC
        
        subplot(2,4,k);
        
        % Plot
        costval = squeeze(METRICS(k,:,i));
        xvalues = MC_fid_inel_xsection / acceptances(k);
        plot(xvalues, costval);
        
        [~,best] = min(costval);
        minsxs(k) = MC_fid_inel_xsection(best);
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
    
    outputfile = sprintf('run_%d_total_fit_%d.pdf', param.runs(1), i);
    print(fig1, sprintf('./combfigs/%s', outputfile), '-dpdf');
    system(sprintf('pdfcrop --margins 10 ./combfigs/%s ./combfigs/%s', outputfile, outputfile));
    
end % Over metrics

end % FIT_ON
