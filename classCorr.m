function icc = classCorr (wNets, inhibCols)

    [trials, rows, cols] = size(wNets);
    e = eye(rows);
    ee = repmat(e, [1 1 trials]);
    eex = ~permute(ee, [3 1 2]);

    
    t = wNets(eex);     %removes diagonal
    t = reshape(t, trials, rows-1, rows);
    tt = t(:, :, (rows-inhibCols+1):rows);
    iNets = reshape(tt, trials, (inhibCols)*(rows-1));
        
        
    [numNets netSize] = size(iNets);
    paircount = 0;
    r = 0;
    
    sigmas = std(iNets');
    means = mean(iNets');
    
    for i = 1:numNets
        
        for j = 1:numNets
        
            if j > i
                
                v1 = iNets(i,:);
                v2 = iNets(j,:);
                m1 = means(i);
                m2 = means(j);
                s1 = sigmas(i);
                s2 = sigmas(j);
                w1 = v1' - m1*ones(length(v1), 1);
                w2 = v2' - m2*ones(length(v2), 1);
                r12 = w1'*w2/(s1*s2*(length(v1) - 1));
                
                
                r = r + r12;
                paircount = paircount +1;
            
            end
        
        end
        
    end
    
    icc = r/paircount;

endfunction
