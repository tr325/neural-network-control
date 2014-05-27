

for i = 1:10
    s = ["~/Documents/Cambridge/NeuralNetwork/src/neural-network-control/stabtest 100 20"];
    system(s);
    W = load("-ascii", "stabilizedW.ascii");
    reEigRec(i,:) = real(eig(W));
end

avgEig = mean(reEigRec);
hist(avgEig)