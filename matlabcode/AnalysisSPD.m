% SPD (silicon pixel detector) fired chip quick quality/outlier analysis
%
% This macro will produce hot and cold chip lists based on kurtosis n-sigma
% cut. Kurtosis as a third-order statistics is very sensitive to outliers.
%
% mikael.mieskolainen@cern.ch, 2018
% Licensed under the MIT License <http://opensource.org/licenses/MIT>.

clear;
close all;

addpath ./src

run    = 274594;
mcname = 'Pythia-6_(17h7a)';

% Detector level
DATA = csvread(sprintf('../figures_xsec/%d/Ascii/TECINT11-B-NOPF-CENTNOTRD_Data-%d_SPDFiredChip.csv', run, run), 0,0);
MC   = csvread(sprintf('../figures_xsec/%d/Ascii/TECINT11-B-NOPF-CENTNOTRD_%s_SPDFiredChip.csv', run, mcname), 0,0);

Nchips = 1200;

DATA = DATA(:,1:Nchips);
MC   = MC(:,1:Nchips);

bits = 0:Nchips-1;


%% Calculate Values Based on Kurtosis

func = @(X) kurtosis(X);

DATA_x = func(DATA);
MC_x   = func(MC);
R      = DATA_x ./ (MC_x);

% Not active
notactive     = isnan(R);
notactive_ind = find(notactive == 1);


%% 2D CHECK THESE OUT, THESE ARE ARBITRARY NOW!!

% Inner Layer
imagesc(reshape(fliplr(DATA_x(1:400)), 40, 10))

% Outer Layer
imagesc(reshape(fliplr(DATA_x(401:1200)), 8, 100))



%%
subplot(3,1,1);
stem(bits, DATA_x,'markersize',0.1)
set(gca,'yscale','log');
title(sprintf('Data Run-%d', run),'interpreter','latex');
ylabel('Kurtosis','interpreter','latex');
axis([-5 1205 0 inf]);

subplot(3,1,2);
stem(bits, MC_x,'markersize',0.1)
set(gca,'yscale','log');
title(sprintf('MC %s', strrep(mcname,'_','-')),'interpreter','latex');
ylabel('Kurtosis','interpreter','latex');
axis([-5 1205 0 inf]);

subplot(3,1,3);
plot(bits, R);
xlabel('SPD Fired Chip','interpreter','latex');
title('DATA/MC','interpreter','latex');
axis([-5 1205 0.5 1.5]);

%%

figure;
hist(R(~notactive), 120);
axis square;
axis([0.5 1.5 0 inf]);
set(gca,'XTick', 0.5:0.1:1.5);
ylabel('Counts','interpreter','latex');
xlabel('Data/MC','interpreter','latex');



%%
% Calculate
medR = median(R(~notactive));
stdR = std(R(~notactive));

fprintf('Median: %0.2f \n', medR);
fprintf('Std:    %0.2f \n', stdR);

xlabel('Data/MC','interpreter','latex');
ylabel('Counts','interpreter','latex');

sigmaCUT = 2.3; % If over this, mark the chip

hotlist  = [];
coldlist = [];
for i = 1:length(R)
    
    % Do not consider completely inactive chips
    if (notactive(i) == 0)
    
        value = R(i) - medR;
        if (abs(value) > sigmaCUT * stdR)
            if (value > 0)
                str = 'HOT';
                hotlist(end+1) = bits(i);
            else
                str = 'COLD';
                coldlist(end+1) = bits(i);
            end
            fprintf('%4.0d : %0.1f sigma \t=> %s \n', bits(i),  value / stdR, str);
        end
    end
end

fprintf('\n\n');

% Print out already MC masked / MC dead
fprintf('Dead/Already masked chips: \n');
fprintf('{');
for i = 1:length(notactive_ind)-1
    fprintf('%d,', notactive_ind(i));
end
fprintf('%d}\n',notactive_ind(end));

% Print out the hot chips
fprintf('\nOutlier chips (HOT): \n');
fprintf('{');
for i = 1:length(hotlist)-1
    fprintf('%d,', hotlist(i));
end
fprintf('%d}\n', hotlist(end));

% Print out the cold chips
fprintf('\nOutlier chips (COLD): \n');
fprintf('{');
for i = 1:length(coldlist)-1
    fprintf('%d,', coldlist(i));
end
fprintf('%d}\n', coldlist(end));

