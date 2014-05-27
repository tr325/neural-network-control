
%% LOADING UP DATA

f1 = 'fEfI8020withplast';
f2 = '.ascii';

%for i = 1:25
%%    
%    kappa0(i, :, :) = load([f1 num2str(i-1) 'Kappa0' f2]);
%    kappa2(i, :, :) = load([f1 num2str(i-1) 'Kappa0.2' f2]);
%    kappa4(i, :, :) = load([f1 num2str(i-1) 'Kappa0.4' f2]);
%    kappa6(i, :, :) = load([f1 num2str(i-1) 'Kappa0.6' f2]);
%    kappa8(i, :, :) = load([f1 num2str(i-1) 'Kappa0.8' f2]);
%    kappa10(i, :, :) = load([f1 num2str(i-1) 'Kappa1' f2]);
%
%end
%refNet = load('fEfI8020withplastREF.asciiKappa');

%%%% eigenspectra
%for i = 1:25
%    
%    k0 = squeeze(kappa0(i, :,:));
%    hold all;figure(1); plot(eig(k0), '.');
%    k2 = squeeze(kappa2(i, :,:));
%    hold all;figure(2); plot(eig(k2), '.');
%    k4 = squeeze(kappa4(i, :,:));
%    hold all;figure(3); plot(eig(k4), '.');
%    k6 = squeeze(kappa6(i, :,:));
%    hold all;figure(4); plot(eig(k6), '.');
%    k8 = squeeze(kappa8(i, :,:));
%    hold all;figure(5); plot(eig(k8), '.');
%    k10 = squeeze(kappa10(i, :,:));
%    hold all;figure(6); plot(eig(k10), '.');
% 
%    %%% max Im and min re parts of eigenspectra
%    mk0(i, 1) = min(real(eig(k0))); mk0(i, 2) = max(imag(eig(k0)));
%    mk2(i, 1) = min(real(eig(k2))); mk2(i, 2) = max(imag(eig(k2)));
%    mk4(i, 1) = min(real(eig(k4))); mk4(i, 2) = max(imag(eig(k4)));
%    mk6(i, 1) = min(real(eig(k6))); mk6(i, 2) = max(imag(eig(k6)));
%    mk8(i, 1) = min(real(eig(k8))); mk8(i, 2) = max(imag(eig(k8)));
%    mk10(i, 1) = min(real(eig(k10))); mk10(i, 2) = max(imag(eig(k10)));
%    
%    
%end

%plotDynamicEvolution(k0);
%plotDynamicEvolution(k2);
%print -depsc2 'kappa2Dyn.eps'
%plotDynamicEvolution(k4);
%plotDynamicEvolution(k6);
%plotDynamicEvolution(k8);
%print -depsc2 'kappa8Dyn.eps'
%plotDynamicEvolution(k10);

%%%%%%%  example eigenspectra plots 
%figure; plot(eig(k0), '.'); set(gca, 'fontsize', 15);
%print -depsc2 'exKappa0eig.eps' 
%figure; plot(eig(k2), '.'); set(gca, 'fontsize', 15);
%print -depsc2 'exKappa2eig.eps' 
%figure; plot(eig(k4), '.'); set(gca, 'fontsize', 15);
%print -depsc2 'exKappa4eig.eps' 
%figure; plot(eig(k6), '.'); set(gca, 'fontsize', 15);
%print -depsc2 'exKappa6eig.eps' 
%figure; plot(eig(k8), '.'); set(gca, 'fontsize', 15);
%print -depsc2 'exKappa8eig.eps' 
%figure; plot(eig(k10), '.'); set(gca, 'fontsize', 15);
%print -depsc2 'exKappa10eig.eps' 


%niceImagesc's of average inhibitory nets
%aK0 = calcAvgInhibNet('fEfI8020withplast', 'Kappa0.ascii', 20, 100);
%aK2 = calcAvgInhibNet('fEfI8020withplast', 'Kappa0.2.ascii', 20, 100);
%aK4 = calcAvgInhibNet('fEfI8020withplast', 'Kappa0.4.ascii', 20, 100);
%aK6 = calcAvgInhibNet('fEfI8020withplast', 'Kappa0.6.ascii', 20, 100);
%aK8 = calcAvgInhibNet('fEfI8020withplast', 'Kappa0.8.ascii', 20, 100);
%aK10 = calcAvgInhibNet('fEfI8020withplast', 'Kappa1.ascii', 20, 100);
%
%niceImagesc(aK0);
%niceImagesc(aK2);
%niceImagesc(aK4);
%niceImagesc(aK6);
%niceImagesc(aK8);
%niceImagesc(aK10);


% ccorr values
kcorr(1) = classCorr(kappa0, 20);
kcorr(2) = classCorr(kappa2, 20);
kcorr(3) = classCorr(kappa4, 20);
kcorr(4) = classCorr(kappa6, 20);
kcorr(5) = classCorr(kappa8, 20);
kcorr(6) = classCorr(kappa10, 20);

lkcorr(1) = classCorr(logical(kappa0), 20);
lkcorr(2) = classCorr(logical(kappa0), 20);
lkcorr(3) = classCorr(logical(kappa4), 20);
lkcorr(4) = classCorr(logical(kappa6), 20);
lkcorr(5) = classCorr(logical(kappa8), 20);
lkcorr(6) = classCorr(logical(kappa10), 20);

disp('kcorr = ')
  kcorr
  
disp('lkcorr = ')
  lkcorr
  
%%%%%%%%  variance calculations %%%%%%%%  
for i = 1:25

    kk0 = squeeze(kappa0(i,:,81:end));
    kk2 = squeeze(kappa2(i,:,81:end));
    kk4 = squeeze(kappa4(i,:,81:end));
    kk6 = squeeze(kappa6(i,:,81:end));
    kk8 = squeeze(kappa8(i,:,81:end));
    kk10 = squeeze(kappa10(i,:,81:end));
    kk0 = kk0(:);
    kk2 = kk2(:);
    kk4 = kk4(:);
    kk6 = kk6(:);
    kk8 = kk8(:);
    kk10 = kk10(:);
    kk0 = kk0(logical(kk0));
    kk2 = kk2(logical(kk2));
    kk4 = kk4(logical(kk4));
    kk6 = kk6(logical(kk6));
    kk8 = kk8(logical(kk8));
    kk10 = kk10(logical(kk10));
    meank0(i) = mean(kk0);
    meank2(i) = mean(kk2);
    meank4(i) = mean(kk4);
    meank6(i) = mean(kk6);
    meank8(i) = mean(kk8);
    meank10(i) = mean(kk10);
    vark0(i) = var(kk0);
    vark2(i) = var(kk2);
    vark4(i) = var(kk4);
    vark6(i) = var(kk6);
    vark8(i) = var(kk8);
    vark10(i) = var(kk10);
end

disp('means')
mean(meank0)
mean(meank2)
mean(meank4)
mean(meank6)
mean(meank8)
mean(meank10)


disp('variances')
mean(vark0)
mean(vark2)
mean(vark4)
mean(vark6)
mean(vark8)
mean(vark10)


%%%%%%%%%  histogrammage %%%%%%%%%% 
%figure; hist(kk0(logical(kk0)), 100); set(gca, 'fontsize', 15);
%print -depsc2 'kappa0hist.eps'
%figure; hist(kk0(kk0(logical(kk0))>-2.5), 100); set(gca, 'fontsize', 15);
%print -depsc2 'kappa0histZoom.eps'
%figure; hist(kk8(logical(kk8)), 100); set(gca, 'fontsize', 15);
%print -depsc2 'kappa8hist.eps'
%figure; hist(kk10(logical(kk10)), 100); set(gca, 'fontsize', 15);
%print -depsc2 'kappa10hist.eps'




%%%%%% reciprocal connection counts %%%%%%%
clk0 = 0;
clk2 = 0;
clk4 = 0;
clk6 = 0;
clk8 = 0;
clk10 = 0;
for i = 1:25

    lk0 = squeeze(logical(kappa0(i, 1:80, 1:80)));
    lk0t = lk0 + lk0';
    clk0 = clk0 + sum(sum(lk0t == 2));
    lk2 = squeeze(logical(kappa2(i, 1:80, 1:80)));
    lk2t = lk2 + lk2';
    clk2 = clk2 + sum(sum(lk2t == 2));    
    lk4 = squeeze(logical(kappa4(i, 1:80, 1:80)));
    lk4t = lk4 + lk4';
    clk4 = clk4 + sum(sum(lk4t == 2));    
    lk6 = squeeze(logical(kappa6(i, 1:80, 1:80)));
    lk6t = lk6 + lk6';
    clk6 = clk6 + sum(sum(lk6t == 2));    
    lk8 = squeeze(logical(kappa8(i, 1:80, 1:80)));
    lk8t = lk8 + lk8';
    clk8 = clk8 + sum(sum(lk8t == 2));
    lk10 = squeeze(logical(kappa10(i, 1:80, 1:80)));
    lk10t = lk10 + lk10';
    clk10 = clk10 + sum(sum(lk10t == 2));

end

disp('Probability of reciprocal connections:')
p0 = clk0/(25*80*80)
p2 = clk2/(25*80*80)
p4 = clk4/(25*80*80)
p6 = clk6/(25*80*80)
p8 = clk8/(25*80*80)
p10 = clk10/(25*80*80)





