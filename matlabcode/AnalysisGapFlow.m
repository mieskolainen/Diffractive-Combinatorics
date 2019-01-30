% GapFlow [EXPERIMENTAL/UNDER WORK]
%
% mikael.mieskolainen@cern.ch, 2018
% Licensed under the MIT License <http://opensource.org/licenses/MIT>.
%
% Combining plots:
% gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=merged.pdf *.pdf

addpath /home/user/cernbox/#Combinatorics_paper/matlab_code
addpath ./src

clear;
close all;

% Number of eigenfunctions
Neig = 5;

% Vector space dimension
N = 6;

% Gapflow type 1,2,3
GPtype = '1';

% Create legend texts
labels = makelegend(6,true);
labels = labels(2:end); % Remove 0-one

% Discretization (CHECK C++ CODE)
MAXVAL = 10.0;

% Run number
run = 274594;

modes = {};
outputfiles = {};
Y = {};

Pmin = 0.0; % Minimum probability

sigma_inel_vis = 72;

% Combination label selection
selected = [];

%%

for source = 1:3
    
    if (source == 1)
        outputfile = sprintf('gapflow_data');
        X = dlmread(sprintf('../figures_xsec/%d/Ascii/gapflow%s_Data-%d.out', run, GPtype, run),' ')';
        labelname = sprintf('Data$_{\\,%d}$', run);
    end
    if (source == 2)
        outputfile = sprintf('gapflow_pythia6');
        X = dlmread(sprintf('../figures_xsec/%d/Ascii/gapflow%s_Pythia-6_(17h7a).out', run, GPtype),' ')';
        labelname = sprintf('Pythia-6$_{\\,%d}$', run);
    end
    if (source == 3)
        outputfile = sprintf('gapflow_phojet');
        X = dlmread(sprintf('../figures_xsec/%d/Ascii/gapflow%s_Phojet_(17h7b).out', run, GPtype),' ')';
        labelname = sprintf('Phojet$_{\\,%d}$', run);
    end
    
    % Create variables
    if (source == 1)
        DISC = size(X,2);
        xval = linspace(0, MAXVAL, DISC);
        entropies = zeros(3,DISC);
    end
    
    modes{source} = labelname;
    outputfiles{source} = outputfile;
    
    % Remove last (empty line)
    if (size(X,1) ~= 2^N)
    X = X(1:end-1,:);
    end
    
    %% Order to Gray code
    
    % Apply GrayCode
    N = 6;
    C_sub = createCBM(N, true);
    %
    graylist = zeros(2^N,0);
    for d = 0:2^N-1
        b = dec2bin(d);
        g = bin2gray(b);
        e = bin2dec(g);
        graylist(d+1) = e; % +1 due to Matlab indexing
        disp({d b g e})
    end
    
    % Order to Graycode
    C_sub_gray = C_sub(graylist'+1,:);
    C_sub = C_sub_gray;
    for j = 1:size(X,2)
        X(:,j) = X(graylist + 1,j); % +1 due to Matlab indexing
    end
    %}
    %figure;
    %stem(X(:,1)); set(gca,'yscale','log');
    
    
    %% Remove first row
    X = X(2:end,:);
    
    %% Plot flow
    %{
    figure;
    plot(COU', '-', 'linewidth', 1.1);
    
    set(gca,'yscale','log');
    axis([0 13 1e-3 1]);
    axis square;
    
    xlabel('$\sum \Delta \eta$','interpreter','latex');
    set(gca,'XTick', round(linspace(0,110,16), 1));
    
    cmd = sprintf('print -dpdf ./gapflowfigs/%s_gapss.pdf', outputfile);
    eval(cmd);
    system(sprintf('pdfcrop --margins 10 ./gapflowfigs/%s_gapss.pdf ./gapflowfigs/%s_gapss.pdf', outputfile, outputfile));
    
    %}
    %
    % Normalize to sum to 1
    %{
    for i = 1:size(X,2)
       X(:,i) = X(:,i) / sum(X(:,i)); 
    end
    %}
    
    Y{source} = X;
    

    %%
    
    
    %cla reset;
    close all;
    figure;
    
    YY = X / sum(X(:,1)) * sigma_inel_vis;
    
    h = plot(xval, X'); hold on;
    %set(gca,'yscale','log');
    %set(gca,'xscale','log');
    
    %axis([0.210 MAXVAL 0 1.10]); axis square;
    xticks(round(linspace(min(xval), max(xval), 11), 1));
    
    xlabel('$z_> = \mathcal{O}_i/\langle \mathcal{O}_i \rangle$','interpreter','latex');
    ylabel('Probability','interpreter','latex');
    title(sprintf('%s: $\\Delta_\\eta$-Combinatorics Flow', labelname),'interpreter','latex');
    
    % Go through legend labels
    %l = legend(labels);
    %set(l, 'interpreter', 'latex', 'fontsize', 10);
    
    if (source == 1)
        for i = 1:length(labels)
            if (X(i, 5) > Pmin)
                selected(end+1) = i;
            end
        end
    end
    
    l = legend(h(selected), labels(selected),'location','northeast'); legend('boxoff');
    set(l,'interpreter','latex','fontsize',7);
    %set(gca,'yscale','log');
    %axis([0 MAXVAL 1e-5 1]);
    %view([90 90]); % Rotate view
    axis square;
    
    cmd = sprintf('print -dpdf ./gapflowfigs/%s.pdf', outputfile);
    eval(cmd);
    system(sprintf('pdfcrop --margins 10 ./gapflowfigs/%s.pdf ./gapflowfigs/%s.pdf', outputfile, outputfile));
    
    
    %% Derivative
    
    close all;
    figure;
    
    YY = X / sum(X(:,1)) * sigma_inel_vis;
    
    D = diff(YY);
    h = plot(xval, D); axis square;
    xlabel('$z_> = \mathcal{O}_i/\langle \mathcal{O}_i \rangle$','interpreter','latex');
    ylabel('$\Delta \sigma_i \equiv \sigma_{i+1} - \sigma_i$ (mb)','interpreter','latex');
    title(sprintf('%s: Hamiltonian component GapFlow', labelname),'interpreter','latex');
    axis([0 MAXVAL -5.0 5.0]);
    
    %view([-90 90]); % Rotate 
    view
    view([90 90]); % Rotate view
    %set(gca,'YDir','reverse');
    set(gca,'YTick', -5:1:5);
    
    %%
    %l = legend(h(selected), labels(selected),'location','northeast'); legend('boxoff');
    %set(l,'interpreter','latex','fontsize',7);
    
    cmd = sprintf('print -dpdf ./gapflowfigs/%s_D.pdf', outputfile);
    eval(cmd);
    system(sprintf('pdfcrop --margins 10 ./gapflowfigs/%s_D.pdf ./gapflowfigs/%s_D.pdf', outputfile, outputfile));
    %}
    
    %{
    figure;
    bar(xval, X','stacked','Barwidth',1,'LineWidth',1e-3); axis square; axis tight;
    
    axis([0 MAXVAL 0 1]); axis square;
    xticks(round(linspace(min(xval), max(xval), 11), 1));
   
    xlabel('$z_> = \mathcal{O}_i/\langle \mathcal{O}_i \rangle$','interpreter','latex');
    
    cmd = sprintf('print -dpng ./gapflowfigs/stacked_%s', outputfile);
    %}
    %eval(cmd);
    %system(sprintf('pdfcrop --margins 10 ./gapflowfigs/stacked_%s ./gapflowfigs/stacked_%s', outputfile, outputfile));
    
    %{
    figure;
    imagesc(xval, 1:21010, log10(X));
    xticks(round(linspace(min(xval), max(xval), 11), 1));
    %}
    % Entropy function
    S = @(p) -sum(p(p>0).*log2(p(p>0)));

    for i = 1:size(X,2)
       entropy(source, i) = S(X(:,i) / sum(X(:,i))); 
    end
    
    % --------------------------------------------------------------------
    % Eigenfunctions
    
    close all;
    
    figure;
    [A,B] = pca(X);
    
    % Create labels
    eig_labels = cell(Neig,1);
    for i = 1:Neig
        eig_labels{i} = sprintf('$e_%d$', i);
    end
    % Plot horizontal line
    h1 = plot(xval([1 end]), zeros(2,1), 'color', ones(3,1)*0.10); hold on;
    
    % Restart color order
    ax = gca;
    ax.ColorOrderIndex = 1;
    
    h2 = plot(xval, A(:,1:Neig), 'linewidth', 1.0); axis square;
    l = legend(h2, eig_labels); legend('boxoff');
    set(l,'interpreter','latex');
    axis([0.210 MAXVAL -0.8 0.8]); axis square;
    xticks(round(linspace(min(xval), max(xval), 11), 1));
    
    xlabel('$z_> = \mathcal{O}_i/\langle \mathcal{O}_i \rangle$','interpreter','latex');
    title(sprintf('%s: $\\Delta_\\eta$-Combinatorics Flow Eigenfunctions', labelname),'interpreter','latex');
    
    cmd = sprintf('print -dpdf ./gapflowfigs/%s_eigen.pdf', outputfile);
    eval(cmd);
    system(sprintf('pdfcrop --margins 10 ./gapflowfigs/%s_eigen.pdf ./gapflowfigs/%s_eigen.pdf', outputfile, outputfile));
   
    
end


%% Factorization

for k = 1:length(Y)

    [F{k}, fact_labels] = factorization(Y{k}, N, labels);
    
    fig = figure;
    plot(xval, F{k}', 'linewidth', 1.0);
    set(gca,'yscale','log');
    axis square;
    axis([1e-5 5 1e-2 1e2]);
    xlabel('$z_> = \mathcal{O}_i/\langle \mathcal{O}_i \rangle$','interpreter','latex');
    xticks(round(linspace(0, max(xval), 11), 1));
    ylabel('$\langle 11 \rangle / (\langle 10\rangle \langle 01 \rangle)^{1/2}$','interpreter','latex');
    title(sprintf('%s: Factorization', modes{k}),'interpreter','latex');
    %l = legend(fact_labels); set(l,'interpreter','latex');
    set(gca,'XTick', 0:0.5:5);
    
    output = sprintf('./gapflowfigs/factorization_%s.pdf', outputfiles{k});
    print(fig, output, '-dpdf');
    system(sprintf('pdfcrop --margins 10 %s %s', output, output));
end


%% Entropy
fig = figure;

plot(xval, entropy, 'linewidth', 1.3); 

xlabel('$z_> = \mathcal{O}_i/\langle \mathcal{O}_i \rangle$','interpreter','latex');
ylabel('Shannon entropy $S$ (bits)','interpreter','latex');

l = legend(modes);
set(l,'interpreter','latex'); legend('boxoff');

axis([0 MAXVAL 0 6]); axis square;
xticks(round(linspace(min(xval), max(xval), 11), 1));
set(gca,'yscale','log');
%set(gca,'xscale','log');

output = sprintf('./gapflowfigs/run_%d_entropy.pdf', run);
print(fig, output, '-dpdf');
system(sprintf('pdfcrop --margins 10 %s %s', output, output));


%% Dispersion

fig = figure;
for index = 1:length(Y)
    
    YY = Y{index} / sum(Y{index}(:,1)) * sigma_inel_vis;
    
    % Dispersion moment
    m = 2;
    YY = diff(YY);
    Hx = (moment(YY, m)).^(1/m);
    
    plot(xval, Hx,'linewidth',1.3); hold on;
    set(gca,'yscale','log');
end
%set(0,'DefaultAxesColorOrder', brewermap(length(labelname),'Dark2'));

xlabel('$z_> = \mathcal{O}_i/\langle \mathcal{O}_i \rangle$','interpreter','latex');
ylabel('Dispersion $\left(\langle \sigma_i^2 \rangle - \langle \sigma_i \rangle^2\right)^{1/2}$ (mb)','interpreter','latex');   
axis square;

%set(gca,'xscale','log');
set(gca,'XTick',[min(xval):1:max(xval)]);
l = legend(modes); set(l,'interpreter','latex');
legend('boxoff');

%axis([0 10 1e-2 10]);
output = sprintf('./gapflowfigs/run_%d_dispersion.pdf', run);
print(fig, output, '-dpdf');
system(sprintf('pdfcrop --margins 10 %s %s', output, output));

