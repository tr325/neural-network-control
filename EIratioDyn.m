% plots dynamic evolution for 50 stabilised W matrices with different 
% numbers of Inhibitory columns

load('stabWs2.mat')

[n rows cols] = size(stabWs);

for i = 1:n 
    W = stabWs(i,:,:);
    W = squeeze(W); 
    if mod(i,5) == 0
      plotDynamicEvolution(W);
      title([num2str(51 - i) ' inhibitory columns']);
    end
    eV = eig(W);
    [d dEInd] = max(imag(eig(W)));
    dEig = eV(dEInd);
    damp(i) = -(real(dEig)-1)/abs(dEig);
    dEigIm(i) = d;
end

%close all
figure
plot(damp)
title('damping factor change w/ EI ratio')

figure
plot(dEigIm, '.')
title('largest Im eignevalue part')