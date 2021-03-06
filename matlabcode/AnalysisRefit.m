% Diffractive Cross Sections Pomeron Intercept and kinematic cutoff refit
%
% mikael.mieskolainen@cern.ch, 2019
% Licensed under the MIT License <http://opensource.org/licenses/MIT>.

clear; close all;
addpath ./src
colors = get(gca, 'ColorOrder');

% **** DISCRETIZATION MANUALLY HERE!! ****
ND = 12;   % Match this with C++ code (one should make this automatic!!)
% **** DISCRETIZATION MANUALLY HERE!! ****

% Ascii file column indices for reading out the fit results
INTERCEPT_ind = 2;
DELTAY_ind    = 3;
COST_ind      = 5;

% #Fit point, Mass re-weight DELTA, Min <DeltaY> cutoff, <neg log(L)>, <KL-divergence>, <KS-error>, <CHI2>
if (COST_ind == 4)
    cost_str = '$-\ln(L)$';
end
if (COST_ind == 5)
    cost_str = 'Kullback-Leibler divergence';
end
if (COST_ind == 6)
    cost_str = 'Kolmogorov-Smirnov error';
end
if (COST_ind == 7)
    cost_str = '$\chi^2$';
end

runs = [274593 274594 274595];
mc   = {'Pythia-6_(17h7a)', 'Phojet_(17h7b)'};

% Create filename
str = '';
for i = 1:length(runs), str = str + "_" + string(runs(i)); end
for i = 1:length(mc),   str = str + "_" + string(mc(i));   end
fp = fopen(sprintf('./fitfigs/refit_%s.tex', str), 'w');

% Create Legends
legs = {};
for i = 1:length(runs), legs{i} = sprintf('Run %d', runs(i)); end


%%

% Over MC models
for model = 1:length(mc)

% Over runs
for run = runs
    
close all;

f1 = figure;
f2 = figure;
figures = {};
i = 1;
    
POM_DELTA = [];

% Over Detector level (1), Unfolded level (2)
for level = 1:2
    
    M = csvread(sprintf('../figures_xsec/%d/Fit/Fit_level_%d_Input_Data-%d_Model_%s.csv', ...
        run, level, run, mc{model}), 1,0);
    
    % Discretization
    x_values = M(1:ND:end, INTERCEPT_ind);
    y_values = M(1:ND, DELTAY_ind);
    X = reshape(M(:, COST_ind), [ND ND]);
    
    % Interpolated values
    FACTOR = 8; % Interpolation factor, one can change it here
    XX = interp2(X, FACTOR, 'spline');
    x_values = linspace(x_values(1),x_values(end), size(XX,2));
    y_values = linspace(y_values(1),y_values(end), size(XX,1));
    
    % Plot
    figures{i} = figure;
    imagesc(x_values, y_values, XX);
    axis square;
    h = colorbar;
    ylabel(h, cost_str,'interpreter','latex');
    
    xlabel('Re-weight $\Delta_P$ (effective Pomeron intercept - 1)','interpreter','latex');
    ylabel('$\langle \Delta Y\rangle_{\min}$ cutoff definition','interpreter','latex');
    set(gca, 'XTick', round(x_values(1):0.01:x_values(end),2));
    set(gca, 'YTick', round(y_values(1):0.5:y_values(end),1));
    title(sprintf('Data-%d, MC: %s', run, strrep(mc{model},'_',' ')), 'interpreter','latex');
    
    % Find optimal point
    min_row = -1;
    min_val = 1e9;
    for k = 1:size(XX,1)
       [a,~] = min(XX(k,:));
       if (a < min_val)
          min_val = a; 
          min_row = k;
       end
    end
    
    % Take slices over x and y-axis
    xx    = XX(min_row,:);
    [a,b] = min(xx);
    yy    = XX(:,b);
    
    POM_DELTA(end+1) = x_values(b);
    fprintf('Model = %s : Run = %d : Effective Pomeron intercept = %0.3f \n', mc{model}, run, POM_DELTA(end));
        
    % Minimum rapidity gap <Delta Y>_min cutoff definition
    figure(f2);
    plot(y_values, yy, 'linewidth', 1.2, 'color', colors(i,:)); axis square; axis tight; hold on;
    set(gca, 'XTick', round(-1e-3:0.5:y_values(end),1) );
    %axis([0 inf min(yy)*0.5 max(yy)*1.05]);
    xlabel('$\langle \Delta Y\rangle_{\min}$ cutoff definition','interpreter','latex');
    ylabel(cost_str,'interpreter','latex');
    
    % Effective Pomeron intercept
    figure(f1);
    plot(x_values, xx, 'linewidth', 1.2, 'color', colors(i,:)); axis square; axis tight; hold on;
    set(gca, 'XTick', round(x_values(1):0.01:x_values(end),2) );
    %axis([-inf inf min(xx)*0.95 max(xx)*1.05]);
    xlabel('Re-weight $\Delta_P$ (effective Pomeron intercept - 1)','interpreter','latex');
    ylabel(cost_str,'interpreter','latex');
    
    % Likelihood mode
    if (COST_ind == 4)
        [min_val, min_ind] = min(xx);
        
        % Now scan the -log(likelihood) + DEVIATION point
        DEVIATION = 0.5; % See formula (6.25) C.Gowan, Statistical Data Analysis
        up_ind = 0;
        for kk = min_ind:1:length(xx)
           if (xx(kk) > min_val + 0.5)
                up_ind = kk;
                break; 
           end
        end
        low_in = 0;
        for kk = min_ind:-1:1
           if (xx(kk) > min_val + 0.5)
                low_in = kk;
                break; 
           end
        end
        
        % Central value
        hold on;
        plot(x_values(min_ind)*ones(5,1),  linspace(min_val, min_val*1.001, 5), 'color', colors(i,:), 'linestyle', '-');
        plot(x_values(low_in)*ones(5,1),   linspace(min_val, min_val*1.001, 5), 'color', colors(i,:), 'linestyle', ':');
        plot(x_values(up_ind)*ones(5,1),   linspace(min_val, min_val*1.001, 5), 'color', colors(i,:), 'linestyle', ':');
        
        %plot(x_values(min_ind)*ones(5,1), linspace( , , 5));
        %text(x_values(min_ind), min_val*1.1, ...
        
        text_str = sprintf('$\\Delta_{P} = %0.3f \\pm %0.3f$', ...
        x_values(min_ind), mean([x_values(min_ind)-x_values(low_in), x_values(up_ind)-x_values(min_ind)]) );
        text(x_values(min_ind)*0.9, ((max(xx) + min(xx))/2)*(1 + 0.0002*i), text_str, 'interpreter','latex', 'fontsize', 15, 'color', colors(i,:));
        
        fprintf('DELTA_P = [%0.5f %0.5f %0.5f] \n', x_values(low_in), x_values(min_ind), x_values(up_ind));
        set(gca,'yscale','log');
    end
    i = i + 1;
    
end % Level

%% Print out

fprintf(fp, '%s & %d & ', mc{model}, run);

for k = 1:length(POM_DELTA)
   fprintf(fp, '%0.3f', POM_DELTA(k));
   if (k < length(POM_DELTA))
       fprintf(fp, ' & ');
   else
       fprintf(fp, ' \\\\ \n');
   end
end
fprintf(fp, '\n');


%% Plot out

for k = 1:length(figures)
figure(figures{k});
outputfile = sprintf('Run_%d_2D_POMERONDELTA_DELTAY_level_%d_model_%d', runs(k), level, model);
print(figures{k}, sprintf('./fitfigs/%s.pdf', outputfile), '-dpdf');
system(sprintf('pdfcrop --margins 10 ./fitfigs/%s.pdf ./fitfigs/%s.pdf', outputfile, outputfile));
end

figure(f1);
l = legend(legs); set(l,'interpreter','latex');
    title(sprintf('MC: %s', strrep(mc{model},'_',' ') ), 'interpreter','latex');
    outputfile = sprintf('POMERONDELTA_level_%d_model_%d', level, model );
    print(f1, sprintf('./fitfigs/%s.pdf', outputfile), '-dpdf');
    system(sprintf('pdfcrop --margins 10 ./fitfigs/%s.pdf ./fitfigs/%s.pdf', outputfile, outputfile));
    
figure(f2);
l = legend(legs); set(l,'interpreter','latex');
    title(sprintf('MC: %s', strrep(mc{model},'_',' ') ), 'interpreter','latex');
    outputfile = sprintf('DELTAY_level_%d_model_%d', level, model );
    print(f2, sprintf('./fitfigs/%s.pdf', outputfile), '-dpdf');
    system(sprintf('pdfcrop --margins 10 ./fitfigs/%s.pdf ./fitfigs/%s.pdf', outputfile, outputfile));

end % runs
end % MC models

fclose(fp);
