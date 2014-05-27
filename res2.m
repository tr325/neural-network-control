
%% LOADING UP DATA
%beta1 = zeros(25, 100, 100);
%beta2 = zeros(25, 100, 100);
%beta3 = zeros(25, 100, 100);
%beta5 = zeros(25, 100, 100);
%beta7 = zeros(25, 100, 100);
f1 = 'fEfI8020withplast';
f2 = '.ascii';

for i = 1:25
    
    beta1(i, :, :) = load([f1 num2str(i-1) 'Beta1' f2]);
    beta2(i, :, :) = load([f1 num2str(i-1) 'Beta2' f2]);

end
for i = 1:4
    
    beta3(i,:,:) = load([f1 num2str(i-1) 'Beta3' f2]);

end

refNet = load('fEfI8020withplastREF.asciiBeta');


%%%% example eigenspectra
%for i = 1:25
%    
%    b1 = squeeze(beta1(i, :,:));
%    hold all;figure(1); plot(eig(b1), '.');
%    b2 = squeeze(beta2(i, :,:));
%    hold all;figure(2); plot(eig(b2), '.');
%    b3 = squeeze(beta3(i, :,:));
%    hold all;figure(3); plot(eig(b3), '.');
%end

for i = 1:25
    
    b1SA(i) = max(real(eig(squeeze(beta1(i,:,:)))));
    b2SA(i) = max(real(eig(squeeze(beta2(i,:,:)))));
    if i < 4
        b3SA(i) = max(real(eig(squeeze(beta3(i,:,:)))));
    end
end

beta2 = beta2(logical(b2SA < 1), :, :);     %beta = 2 set
%beta3 = beta3(1:7, :, :);       %the simulation broke after this


%Example eigenvalue plots
figure; hold on; 
plot(eig(refNet), 'r.');
plot(eig(squeeze(beta1(1,:,:))), '.');
axis('equal'); axis('square');
set(gca, 'fontsize', 15)
print -depsc2 'beta1eigs.eps'

figure; hold on;
plot(eig(2.*refNet), 'r.');
plot(eig(squeeze(beta2(1,:,:))), '.');
axis('equal'); axis('square');
set(gca, 'fontsize', 15)
print -depsc2 'beta2eigs.eps'

figure; hold on;
plot(eig(3.*refNet), 'r.');
plot(eig(squeeze(beta3(1,:,:))), '.');
axis('equal'); axis('square');
set(gca, 'fontsize', 15)
print -depsc2 'beta3eigs.eps'

%some dynamics plots
%plotDynamicEvolution(squeeze(beta1(1,:,:)));
%print -depsc2 'beta1dyn.eps'
%plotDynamicEvolution(squeeze(beta2(1,:,:)));
%print -depsc2 'beta2dyn.eps'
%plotDynamicEvolution(squeeze(beta3(1,:,:)));
%print -depsc2 'beta3dyn.eps'
%
%%Nice imagescs of example networks
%niceImagesc(squeeze(beta1(1,:,:)), 20);
%print -depsc2 'exBeta1net.eps'
%niceImagesc(squeeze(beta2(1,:,:)), 20);
%print -depsc2 'exBeta2net.eps'
%niceImagesc(squeeze(beta3(1,:,:)), 20);

bcorr(1) = classCorr(beta1, 20);
bcorr(2) = classCorr(beta2, 20);
bcorr(3) = classCorr(beta3, 20); 

lbcorr(1) = classCorr(logical(beta1), 20);
lbcorr(2) = classCorr(logical(beta2), 20);
lbcorr(3) = classCorr(logical(beta3), 20);

%Nice imagescs of agerage nets. 
ab1 = calcAvgInhibNet('fEfI8020withplast', 'Beta1.ascii', 20, 100);
ab2 = calcAvgInhibNet('fEfI8020withplast', 'Beta2.ascii', 20, 100);
niceImagesc(ab1);
print -depsc2 'avgB1net.eps'
niceImagesc(ab2);
print -depsc2 'avgB2net.eps'


