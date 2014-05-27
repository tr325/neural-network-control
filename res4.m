f1 = 'fEfI8020withplast';
f2 = '.ascii';

for i = 1:25

    kappa(i, :, :) = load([f1 num2str(i-1) 'Kappa0.4' f2]);

end

refNet = load('exKappaRefNet.ascii');

% example eigenspectrum
k1 = squeeze(kappa(2,:,:));
figure; plot(eig(k1), '.'); 
set(gca, 'fontsize', 13)
print -depsc2 'exKappaEigs.eps'

% some dynamics. woopwoop. 
plotDynamicEvolution(k1);
print -depsc2 'exKappaDyn.eps'


% imgesc heatmaps
aK4 = calcAvgInhibNet('fEfI8020withplast', 'Kappa0.4.ascii', 20, 100);
niceImagesc(aK4);
print -depsc2 'avgKappa0.4Net.eps'
niceImagesc(k1(:,81:end));
print -depsc2 'exKappa0.4Net.eps'


% corrrrrelation values
icc = classCorr(kappa, 20)
lcc = classCorr(logical(kappa), 20)










