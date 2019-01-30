% Unfolding Regularization (number of EM-iterations) Systematics
%
% mikael.mieskolainen@cern.ch, 2018
% Licensed under the MIT License <http://opensource.org/licenses/MIT>.

close all; clear;
addpath ./src

for run = [274593]% 274594 274595]

    fig = {};
    for i = 1:4
       fig{i} = figure; 
    end    
    modelname = {'Pythia 6','Phojet'};
    
    for model = 1:2
    
    % Vector space dimension
    N = 6;

    % Detector level
    X_raw = csvread(sprintf('../figures_xsec/%d/Ascii/Data-%d_x_rates.csv', run, run), 1,0);
    X_raw = X_raw(2:end,2); % Take only fiducial

    % Unfolded
    iters = 1:25;
    X = zeros(2^N-1,length(iters));
    k = 1;
    for iter = iters
        M = csvread(sprintf('../figures_xsec/%d/Ascii/Data-%d_x_unfolded_rates_iter_%d_model_%d.csv', run, run, iter, model), 1,0);
        X(:,k) = M(2:end,2); % Take only fiducial (not extrapolated / 0-bin)
        k = k + 1;
    end

    % Normalize all
    X_raw = X_raw / sum(X_raw);
    
    for i = 1:size(X,2)
       X(:,i) = X(:,i) / sum(X(:,i));
    end

    % Calculate entropies
    H_raw = -sum( X_raw(X_raw > 0) .* log2(X_raw(X_raw > 0)));
    H = zeros(1,size(X,2));
    for i = 1:size(X,2)
       H(i) = -sum( X(:,i) .* log2(X(:,i))); 
    end

    % Calculate relative entropy
    for i = 1:size(X,2)
       KL(i) = sum( log2(X_raw ./ (X(:,i) + 1e-12) + 1e-12) .* X_raw ); 
    end

    %% Calculate covariance matrix

    % Normalize to cross sections
    for i = 1:size(X,2)
       X(:,i) = X(:,i) / sum(X(:,i)) * 72;
    end

    CC = cov(X');
    figure;
    imagesc( CC ); colorbar;
    axis square;

    %% Print out relative 1 sigma uncertainty

    relerr = diag(CC) ./ X(:,5) * 100;
    [(1:63)', relerr]


    %%
    figure(fig{1});
    plot(KL, 's-'); hold on;
    set(gca,'XTick', [0 iters]); axis square;
    xlabel('Regularization (unfolding iterations)','interpreter','latex');
    ylabel('Relative entropy $KL( P_{raw}(\Omega_k) | P_{unfold}(\Omega_k)) $ (bits)','interpreter','latex');

    % Plot detector level and unfolded
    figure(fig{2});
    if (model == 1) % only once
    plot(0, H_raw, 'ro'); hold on; % detector level
    end
    plot(iters, H, 's-'); hold on;
    set(gca,'XTick', [0 iters]); axis square;
    xlabel('Regularization (unfolding iterations)','interpreter','latex');
    ylabel('$H( P(\Omega_k) $ (bits)','interpreter','latex');
    %set(gca,'yscale','log');
    %set(gca,'xscale','log');
    
    % Plot
    figure(fig{3});
    plot(KL, H, 's-'); axis square; hold on;
    for i = 1:length(H)
        text(KL(i), H(i), sprintf('%d', i), 'interpreter','latex', 'fontsize', 12);
    end
    xlabel('Relative entropy $KL( P_{raw}(\Omega_k) | P_{unfold}(\Omega_k)) $ (bits)','interpreter','latex');
    ylabel('Entropy $H( P_{unfold}(\Omega_k)) $ (bits)','interpreter','latex');
    set(gca,'yscale','log');
    %set(gca,'xscale','log');
    
    % Minimum Description Length
    figure(fig{4});
    plot(1:25, H + KL, 'o-'); axis square; hold on;
    ylabel('$H( P(\Omega_k) + KL( P_{raw}(\Omega_k) | P_(\Omega_k)) $ (bits)','interpreter','latex');
    xlabel('Regularization (unfolding iterations)','interpreter','latex');
    set(gca,'XTick', [0 iters]);
    %set(gca,'xscale','log');
    
    % Plot the running of all combinations
    figX = figure;
    plot(X'); set(gca,'yscale','log'); axis square;
    set(gca,'XTick', [0 iters]);
    xlabel('Regularization (unfolding iterations)',  'interpreter', 'latex');
    ylabel('Partial cross section $\sigma_k$ (mb)',  'interpreter', 'latex');
    title(sprintf('Run: %d, Model: %s', run, modelname{model}), 'interpreter', 'latex');
    axis([0 25 1e-4 1e2]);
    
    outputfile = sprintf('Run_%d_unfold_running_model_%d', run, model);
    print(figX, sprintf('./unfoldfigs/%s.pdf', outputfile), '-dpdf');
    system(sprintf('pdfcrop --margins 10 ./unfoldfigs/%s.pdf ./unfoldfigs/%s.pdf', outputfile, outputfile));
    
    
    end % Model loop
    
    for i = 1:length(fig)
        figure(fig{i});
        l = legend(modelname);
        set(l,'interpreter','latex');
    end
    
    outputfile = sprintf('Run_%d_unfold_iterations', run);
    print(fig{1}, sprintf('./unfoldfigs/%s.pdf', outputfile), '-dpdf');
    system(sprintf('pdfcrop --margins 10 ./unfoldfigs/%s.pdf ./unfoldfigs/%s.pdf', outputfile, outputfile));
    
    outputfile = sprintf('Run_%d_unfold', run);
    print(fig{2}, sprintf('./unfoldfigs/%s.pdf', outputfile), '-dpdf');
    system(sprintf('pdfcrop --margins 10 ./unfoldfigs/%s.pdf ./unfoldfigs/%s.pdf', outputfile, outputfile));
    
    outputfile = sprintf('Run_%d_unfold_entropy', run);
    print(fig{3}, sprintf('./unfoldfigs/%s.pdf', outputfile), '-dpdf');
    system(sprintf('pdfcrop --margins 10 ./unfoldfigs/%s.pdf ./unfoldfigs/%s.pdf', outputfile, outputfile));
    
    outputfile = sprintf('Run_%d_unfold_MDL', run);
    print(fig{4}, sprintf('./unfoldfigs/%s.pdf', outputfile), '-dpdf');
    system(sprintf('pdfcrop --margins 10 ./unfoldfigs/%s.pdf ./unfoldfigs/%s.pdf', outputfile, outputfile));
    
    
end % Run loop



