% processes results of scaletest from c++ for visiualisations

Wref = load("-ascii", "Wref.ascii");
W1 = load("-ascii", "Wbeta1.ascii");
W2 = load("-ascii", "Wbeta2.ascii");
W4 = load("-ascii", "Wbeta4.ascii");
W8 = load("-ascii", "Wbeta8.ascii");

inhibCols = 20;

I1 = W1(:, 80:end);
I2 = W2(:, 80:end);
I4 = W4(:, 80:end);
I8 = W8(:, 80:end);

I1 = I1(:);
I2 = I2(:);
I4 = I4(:);
I8 = I8(:);

in1 = abs(I1(find(I1)));
in2 = abs(I2(find(I2)));
in4 = abs(I4(find(I4)));
in8 = abs(I8(find(I8)));

a = length(in1);
if length(in2) < a
    a = length(in2);
end
if length(in4) < a
    a = length(in4);
end
if length(in8) < a
    a = length(in8);
end

in1short = in1(1:a);
in2short = in2(1:a);
in4short = in4(1:a);
in8short = in8(1:a);

figure
hold all
plot(in1short, in4short, '.');
plot(in2short, in4short, '.');
plot(in4short, in4short, '.');
plot(in8short, in4short, '.');
ylabel('Inhibitory weights for Beta = 4', 'fontsize', 15);
xlabel('Inhibitory weght values for other values of beta', 'fontsize', 15);
legend('B = 1', 'B = 2', 'B = 4', 'B = 8')

