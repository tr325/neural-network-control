
## Copyright (C) 2014 Tom Renner

%plots a load of SSA, SA, and stuff against E:I ratio

    dim = 100;
    eigRadius = [];
    smoothSAStart = [];
    smoothSAEnd = [];
    saStart = [];
    saEnd = [];
    stabWs = [];        %stores the stabilised W matrices

    
    for i = 0:40
        in = 50-i;
        s = ["~/Documents/Cambridge/NeuralNetwork/src/neural-network-control/stabtest 100 " num2str(in)];
        system(s);
        W1 = load("-ascii", "generatedW.ascii");
        W2 = load("-ascii", "stabilizedW.ascii");
        saStart = [saStart; max(real(eig(W1)))];
        saEnd = [saEnd; max(real(eig(W2)))];
        eigRadius = [eigRadius; max(abs(eig(W1)))];
        SSA = load("-ascii", "SSArecord.ascii");
        smoothSAStart = [smoothSAStart; SSA(1, 1)];
        smoothSAEnd = [smoothSAEnd; SSA(1, 2)];
        stabWs(i+1, :, :) = W2;
    end
    
    DATA2  = [eigRadius saStart saEnd smoothSAStart smoothSAEnd];
    save("-ascii", "SSAvsEIdata3.ascii", "DATA");
    
    save('stabWs2.mat','stabWs')
    figure
    hold all
    for i = 1:5
        plot(DATA(:, i), '.')
    end
