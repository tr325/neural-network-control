function [IAvg EAvg] = IEstatistics(W, inhibCols)
    
    [rows SIZE] = size(W);
    exCols = SIZE - inhibCols;
    
    eW = W(:, 1:exCols);
    iW = W(:, exCols:end);
    
    [eIdx j eWstrengths] = find(eW(:));
    [iIdx j iWstrengths] = find(iW(:));
    
    EAvg = sum(eWstrengths)/length(eWstrengths);
    IAvg = sum(iWstrengths)/length(iWstrengths);
    
end
  