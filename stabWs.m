%runs stabtest many times and stores the stabilised matrices

stabWs = [];

for i = 0:40
    in = 50-i;
    s = ["~/Documents/Cambridge/NeuralNetwork/src/neural-network-control/stabtest 100 " num2str(in)];
    system(s);
    sW = load("-ascii", "stabilizedW.ascii");
    stabWs(i) = sW;
end