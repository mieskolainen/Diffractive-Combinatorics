% Plot folding matrices
%
% mikael.mieskolainen@cern.ch, 2018
% Licensed under the MIT License <http://opensource.org/licenses/MIT>.

clear; close all;
addpath /home/user/cernbox/#matlabcodes
addpath ./src

basepath = '/home/user/cernbox/ALICE/Diffractive-Combinatorics/figures_xsec';
run = 274595;

% Monte Carlos
mctypes = {'Pythia-6_(17h7a)', 'Phojet_(17h7b)'};

for k = 1:length(mctypes)

fig1 = figure('units', 'normalized', 'outerposition', [0 0 1 1]);

mc = mctypes{k};

filename = sprintf('%s/%d/Unfolding/FoldingMatrix_%s.csv', basepath, run, mc);
X = csvread(filename);

% Normalize columns to 1
for i = 1:size(X,2)
   X(:,i) = X(:,i) / sum(X(:,i)); 
end
X = X'; % Transpose
N = log2(size(X,1));

X = X * 100; % To percentage

imagesc(0:2^N-1, 0:2^N-1, X); hold on;

THRESH = 1;
for i = 1:size(X,1)
    for j = 1:size(X,2)
        
        value = X(j,i);
        if (value > THRESH)
        space = '';
        if (round(value) < 10)
            space = ' ';
        end
        
        color = [0 0 0];
        if (value > 50)
        color = [1 1 1];
        end
        
        text(i-1.4, j-1, sprintf('%s%0.0f', space,value),'color',color,'fontsize', 5);
        end
    end
end

%set(gca,'YDir','normal');
set(gca,'FontSize', 6);

xlabel('Generator', 'interpreter', 'latex', 'fontsize', 12);
ylabel('Detector',  'interpreter', 'latex', 'fontsize', 12);

%labels = {};
%for i = 0:size(X,1)-1
%   labels{end+1} = sprintf('%d', i); 
%end
labels = makelegend(N, true);

xticks(0:size(X,1)-1);
yticks(0:size(X,1)-1);
xticklabels(labels); xtickangle(90);
yticklabels(labels);
set(gca,'TickLabelInterpreter', 'latex');

h = gca;
h.XAxis.TickLength = [0 0];
h.YAxis.TickLength = [0 0];

axis square;
colormap(flipud(hot(256)));
%colormap(hot(256));
%colorbar;
caxis([0 100]);

outputfile = sprintf('Run_%d_foldingmatrix_MC_%d', run, k);
print(fig1, sprintf('./unfoldfigs/%s.pdf', outputfile), '-dpdf');
system(sprintf('pdfcrop --margins 2 ./unfoldfigs/%s.pdf ./unfoldfigs/%s.pdf', outputfile, outputfile));

end

%%
system('source copyplots.sh');

