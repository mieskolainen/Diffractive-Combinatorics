%{
TriggerData::Print:: .Name = CINT11-B-NOPF-CENTNOTRD, .BCMask = B, .LMa = 464881908.00, .LMb = 0.00, .L0aL0b = 0.66
TriggerData::Print:: .Name = CINT11-I-NOPF-CENTNOTRD, .BCMask = I, .LMa = 46488191.00, .LMb = 0.00, .L0aL0b = 0.67
TriggerData::Print:: .Name = CINT11-A-NOPF-CENTNOTRD, .BCMask = A, .LMa = 860031532.00, .LMb = 0.00, .L0aL0b = 0.17
TriggerData::Print:: .Name = CINT11-C-NOPF-CENTNOTRD, .BCMask = C, .LMa = 836787440.00, .LMb = 0.00, .L0aL0b = 0.19
TriggerData::Print:: .Name = CINT11-E-NOPF-CENTNOTRD, .BCMask = E, .LMa = 46488192.00, .LMb = 0.00, .L0aL0b = 0.50
%}

(464881908.00 / 860031532.00) * (0.66 / 0.17)

%% Plot

%{
close all;
errorbar(1:2^N-1, rates / sum(rates) * sigma_inel_fiducial, sqrt(rates) / sum(rates) * sigma_inel_fiducial, 'k.', 'markersize', 6);
hold on;

% Count rates
for i = 1:2
   x = zeros(2^N-1,1);
   for c = 1:2^N-1
      x(c) = sum(X{i}(:,comb_level) == c);
   end
   errorbar(1:2^N-1, x / sum(x) * sigma_inel_fiducial, sqrt(x) / sum(x) * sigma_inel_fiducial, '.', 'markersize', 6);
end
ylabel('mb','interpreter','latex');
set(gca,'yscale','log');
axis tight;
set(gca,'XTick',1:2^N-1);
%}
