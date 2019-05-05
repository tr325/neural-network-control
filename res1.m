% first set of report results

fNames = cell(4);
fnames = {'fixedE5050noplast', 'fixedE8020noplast', 'fixedE5050withPlast', 'fixedE8020withPlast'};

%
resWs = zeros(100, 100, 100);
%
for n = 1:4
    for i = 1:25
        
        f = [fnames{n}, num2str(i-1), '.ascii'];
        W = load(f);
        resWs((n-1)*25 + i, :, :) = W;
        
    end
end
%
nop5050s = resWs(1:25, :, :);
nop8020s = resWs(26:50, :, :);
wp5050s = resWs(51:75, :, :);
wp8020s = resWs(76:100, :, :);

nop5050Is = nop5050s(:, :, 51:end);
nop8020Is = nop8020s(:, :, 81:end);
wp5050Is = wp5050s(:, :, 51:end);
wp8020Is = wp8020s(:, :, 81:end);
%
resIs = resWs(:, :, 81:end);
%

ccorr(1) = classCorr(nop5050s, 50);
ccorr(2) = classCorr(nop8020s, 20);
ccorr(3) = classCorr(wp5050s, 50);
ccorr(4) = classCorr(wp8020s, 20);



    
disp('Correlation values for 5050np, 8020np, 5050wp, 8020wp:')
  ccorr


% correlation values for synapse locations (only, not strengths) of each set. 
lnp5050s = logical(nop5050s);
lnp8020s = logical(nop8020s);
lwp5050s = logical(wp5050s);
lwp8020s = logical(wp8020s);

lcorr(1) = classCorr(lnp5050s, 50);
lcorr(2) = classCorr(lnp8020s, 20);
lcorr(3) = classCorr(lwp5050s, 50);
lcorr(4) = classCorr(lwp8020s, 20);

disp('Correlation values for LOGICAL 5050np, 8020np, 5050wp, 8020wp:')
  lcorr
 
%Making average network images
a5050np = calcAvgInhibNet(fnames{1}, '.ascii', 50, 100);
%niceImagesc(a5050np);
%%print -depsc2 'avgfE5050np.eps'
a8020np = calcAvgInhibNet(fnames{2}, '.ascii',20, 100);
%niceImagesc(a8020np);  
%%print -depsc2 'avgfE8020np.eps'
a5050wp = calcAvgInhibNet(fnames{3},'.ascii', 50, 100);
%niceImagesc(a5050wp); 
%%print -depsc2 'avgfE5050wp.eps'
a8020wp = calcAvgInhibNet(fnames{4}, '.ascii',20, 100);
%niceImagesc(a8020wp);
%print -depsc2 'avgfE8020wp.eps'


%eigenvalue plots - GET A REFERENCE TO COMPARE THEM TO
%p1 = squeeze(nop5050s(15,:,:));
%figure; plot(eig(p1), '.');
%axis('equal'); axis('square');
%%title('np5050')
%print -depsc2 'fE5050npEigs.eps'
%p2 = squeeze(nop8020s(12,:,:));
%figure; plot(eig(p2), '.');
%%title('np8020')
%axis('equal'); axis('square');
%print -depsc2 'fE8020npEigs.eps'
%p3 = squeeze(wp5050s(12,:,:));
%figure; plot(eig(p3), '.');
%%title('wp5050')
%axis('equal'); axis('square');
%print -depsc2 'fE5050wpEigs.eps'
%p4 = squeeze(wp8020s(12,:,:));
%figure; plot(eig(p4), '.');
%%title('wp8020')
%axis('equal'); axis('square');
%print -depsc2 'fE8020wpEigs.eps'


% plots of dynamics of these matrices
%plotDynamicEvolution(p1)
%print -depsc2 'fE5050npDyn.eps'
%plotDynamicEvolution(p2)
%print -depsc2 'fE8020npDyn.eps'
%plotDynamicEvolution(p3)
%print -depsc2 'fE5050wpDyn.eps'
%plotDynamicEvolution(p4)
%print -depsc2 'fE8020wpDyn.eps'


