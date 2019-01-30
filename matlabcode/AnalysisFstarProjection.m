% F*-projection algorithm for pseudorapidity gap distributions
%
% mikael.mieskolainen@cern.ch, 2018
% Licensed under the MIT License <http://opensource.org/licenses/MIT>.

clear; close all;
addpath /home/user/cernbox/#matlabcodes
addpath ./src

basepath = '/home/user/cernbox/ALICE/offline_new/figures_xsec';

% 274593, 274594, 274595
run = 274595;

% **** THESE COULD BE MORE AUTOMATIC, NOW MANUAL!! *****

N = 6; % Vector space dimension

% Fiducial unfolded inelastic in DATA (mb)
sigma_inel_fiducial = 71.5;

% Generator level fiducial in MC generators (mb)
sigma_inel_fiducial_MC = [76 77];

% Minimum and Maximum eta
etamin = -7.0;
etamax =  6.3;

%comb_level = 1;  % Detector level
comb_level = 2;   % Generator (Fiducial) level

% This is raw rates
%DATA = csvread(sprintf('%s/%d/Ascii/TECINT11-B-NOPF-CENTNOTRD_Data-%d.csv', basepath, run, run), 1,0);

% This is detector level
if (comb_level == 1)
    DATA = csvread(sprintf('%s/%d/Ascii/Data-%d_x_rates.csv', basepath, run, run), 1,0);
end
% This is beam-gas substracted and unfolded rates -> generator (fiducial) level
if (comb_level == 2)
    unfold_model = 1; % Unfolding model
    unfold_iter  = 5; % Number of EM-iterations
    DATA = csvread(sprintf('%s/%d/Ascii/Data-%d_x_unfolded_rates_iter_%d_model_%d.csv', basepath, run, run, unfold_iter, unfold_model), 1,0);
end

% MC generators
X{1} = csvread(sprintf('%s/%d/Ascii/TECINT11-B-NOPF-CENTNOTRD_Pythia-6_(17h7a).csv', basepath, run), 1,0);
X{2} = csvread(sprintf('%s/%d/Ascii/TECINT11-B-NOPF-CENTNOTRD_Phojet_(17h7b).csv',   basepath, run), 1,0);

dataname  = sprintf('Data$_{\\,%d}$ $F^{*}$-projected', run);
mcname{1} = 'Pythia 6';
mcname{2} = 'Phojet';

for ANALYSIS = 1:5

fprintf('Analysis %d/%d \n', ANALYSIS, 5);

% Variable index
if (ANALYSIS == 1) % Negative side gap (backward)
    variable_index = 10;
    eta_edge = -7; % ADC left edge
end
if (ANALYSIS == 2) % Positive side gap (forward)
    variable_index = 11;
    eta_edge = 6.3; % ADA right edge
end
if (ANALYSIS == 3) % Largest gap of the event
    variable_index = 12;
    eta_edge = 0;
end

%% Construct abstract frame basis based on MC

binedges  = linspace(0, abs(etamax - etamin), 20);
bincenter = (binedges(2:end) + binedges(1:end-1))/2;
DELTA     = (binedges(2) - binedges(1));

densities = cell(2,1);
for i = 1:2
    densities{i} = zeros(2^N-1,length(binedges)-1);
    
    % Loop over combinations
    for c = 1:2^N-1
        
        if (ANALYSIS == 1 || ANALYSIS == 2 || ANALYSIS == 3)
            gap = abs( X{i}(X{i}(:,comb_level) == c, variable_index) - eta_edge );
        end
        if (ANALYSIS == 4)
            gap = max(abs( X{i}(X{i}(:,comb_level) == c, 10) - etamin), abs( X{i}(X{i}(:,comb_level) == c, 11) - etamax )  );
        end
        if (ANALYSIS == 5)
            gap = max(max(abs( X{i}(X{i}(:,comb_level) == c, 10) - etamin), abs( X{i}(X{i}(:,comb_level) == c, 11) - etamax )), X{i}(X{i}(:,comb_level) == c, 12) );
        end
        
        a = hist1m(gap, binedges);
        a = a / sum(a);
        densities{i}(c,:) = a;
        %stephistedge(binedges,a); hold on;
        
        gap1 = abs( X{i}(X{i}(:,comb_level) == c, 10) - etamin );
        gap2 = abs( X{i}(X{i}(:,comb_level) == c, 11) - etamax );
        
        CC = hist2m([gap1 gap2], binedges, binedges);
    end
end


%% REAL DATA

fig1 = figure;

% Count 2^N-1 rates
rates = DATA(2:end, 2);

% Construct rapidity gap distribution

% First plot MC
for i = 1:2
    
    ind = (X{i}(:,comb_level) ~= 0); % Choose non-zero

    if (ANALYSIS == 1 || ANALYSIS == 2 || ANALYSIS == 3)
        gap = abs( X{i}(ind, variable_index) - eta_edge );
    end
    if (ANALYSIS == 4)
        gap = max(abs( X{i}(ind, 10) - etamin), abs( X{i}(ind, 11) - etamax ) );
    end
    if (ANALYSIS == 5)
        gap = max( max(abs( X{i}(ind, 10) - etamin), abs( X{i}(ind, 11) - etamax )), X{i}(ind, variable_index) );
    end
    
    a = hist1m(gap, binedges);
    a = a / sum(a);
    stephistedge(binedges, sigma_inel_fiducial_MC(i) * a / DELTA, 'linewidth', 1); hold on;
end

% DATA
gapdensity = cell(2,1);
for i = 1:2
    gapdensity{i} = zeros(1, length(binedges)-1);
    for c = 1:2^N-1
       gapdensity{i} = gapdensity{i} + rates(c)*densities{i}(c,:); 
    end
end

% Take the mean of different different model projections
x = (gapdensity{1}+gapdensity{2})/2;
x_err = sqrt(x);

% Plot
errorbar(bincenter, x * sigma_inel_fiducial / DELTA / sum(x), ...
    x_err * sigma_inel_fiducial / DELTA / sum(x), ...
    x_err * sigma_inel_fiducial / DELTA / sum(x), 'ks', 'markersize', 1.5, 'CapSize', 1); hold on;

% Combined Systematic, which takes into account all uncertainties
syst_up = 1.02*max([gapdensity{1}; gapdensity{2}] / sum(rates));
syst_do = 0.98*min([gapdensity{1}; gapdensity{2}] / sum(rates));
stepfilledge(binedges, syst_up * sigma_inel_fiducial / DELTA, ...
                       syst_do * sigma_inel_fiducial / DELTA, [0 0 0],[0 0 0],0.1);

if (ANALYSIS == 1)
    xlabel(sprintf('$\\Delta\\eta^B \\equiv \\eta_{min} + %0.1f$', -etamin), 'interpreter','latex');
    ylabel('$\frac{d\sigma}{d\Delta\eta^{B}}$ (mb)','interpreter','latex');
end
if (ANALYSIS == 2)
    xlabel(sprintf('$\\Delta\\eta^F \\equiv |\\eta_{max} - %0.1f|$', etamax), 'interpreter','latex');
    ylabel('$\frac{d\sigma}{d\Delta\eta^{F}}$ (mb)', 'interpreter','latex');
end
if (ANALYSIS == 3)
    xlabel('$\Delta\eta^{D} \equiv$ max($|\eta_{i} - \eta_{i+1}|$)','interpreter','latex');
    ylabel('$\frac{d\sigma}{d\Delta\eta^{D}}$ (mb)', 'interpreter','latex');
end
if (ANALYSIS == 4)
    xlabel('max(\{$\Delta\eta^{B}, \Delta\eta^{F}$\})','interpreter','latex');
    ylabel('$\frac{d\sigma}{d\Delta\eta}$ (mb)', 'interpreter','latex');
end
if (ANALYSIS == 5)
    xlabel('max(\{$\Delta\eta^{B}, \Delta\eta^{F}, \Delta\eta^{D}$\})','interpreter','latex');
    ylabel('$\frac{d\sigma}{d\Delta\eta}$ (mb)', 'interpreter','latex');
end

title(sprintf('$%0.1f < \\eta < %0.1f, \\, p_t > 0.05$ GeV, charged particle level', etamin, etamax), 'interpreter', 'latex');

axis square; set(gca,'yscale','log');
l = legend([mcname, dataname]); set(l,'interpreter','latex','fontsize',12);legend('boxoff');

set(gca, 'XTick', 0:1:13);

if (ANALYSIS ~= 3)
axis([0 13.0 1e-2 3e2]);
else
axis([0 13.0 1e-3 3e2]);
end
outputfile = sprintf('Run_%d_gaptype_%d', run, ANALYSIS);
print(fig1, sprintf('./gapfigs/%s.pdf', outputfile), '-dpdf');
system(sprintf('pdfcrop --margins 10 ./gapfigs/%s.pdf ./gapfigs/%s.pdf', outputfile, outputfile));


end

