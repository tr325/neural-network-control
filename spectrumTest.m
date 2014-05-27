

dim = 100; 
eigRadius = [];

for i = 1:50
    in = 50-i+1;
    s = ["~/Documents/Cambridge/NeuralNetwork/src/neural-network-control/gentest 100 " num2str(in)];
    system(s);
    genW = load("-ascii", "generatedTestW.ascii");
    genB = load("-ascii", "generatedTestB.ascii");
    
%    figure
%    imagesc(genB)
%    colorbar

    eigRadius = [eigRadius max(abs(eig(genW)))];
end

figure
plot(eigRadius)
