function [avgINet] = calcAvgInhibNet (fname1, fname2, inhibCols, SIZE)
    %Calculates the average inhibitory network of the set of 25 data files with the prefix "fname"
    
    x = SIZE - inhibCols;   %index of first inhibitory column 
    sumINets = zeros(SIZE, SIZE - x);  %container for sum of all inhibitory nets
    
    for i = 0:24
        
        f = [fname1, num2str(i), fname2];
        W = load(f, '-ascii');
        W(logical(eye(length(W)))) = 0;
        sumINets = sumINets + W(:, x+1:end);
    
    end
  
    avgINet = sumINets/25;        
        
end
