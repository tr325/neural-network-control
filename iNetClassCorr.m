function [icc] = iNetClassCorr (fname, inhibCols, SIZE)
    % Calculates the intra-class correlation for a set of simulation results 
    % ie. the generalisation of the correlation coefficient of a pair of 
    % frames to that for a set of many frames. 
    
    x = SIZE - inhibCols;   %index of first inhibitory column 
    iNets = [];
    
    for i = 0:24
        
        f = [fname, num2str(i), ".ascii"];
        W = load(f, '-ascii');
        I = W(:, x+1:end);
        iNets = [iNets, I(:)];
        
    end   
    
    icc = classCorr(iNets);
    
    
    figure;
    imagesc(W(:, x+1:end));
    colorbar
    title([fname ' sample inhibition network']);
    m = mean(iNets, 2);
    avgNet = reshape(m, SIZE, inhibCols);
    figure;
    imagesc(avgNet);
    colorbar
    title([fname ' average inhibition network']);

    
    
    
%    [K N] = size(iNets)
%    xBar = sum(sum(iNets))/(K*N)
%    s2 = 1/(N*K)*(sum(sum((iNets - xBar*ones(K, N)).^2)))
%    xNBars = mean(iNets)
%    
%    icc = sum((xNBars - xBar*ones(1, N)).^2)/(N*s2)
    
    
endfunction
