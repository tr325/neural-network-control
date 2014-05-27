function niceImagesc (net)

    [nr nc] = size(net);
    net(net == 0) = nan;    %pcolor can be set to ignore nans
    
    figure;
    h = pcolor([flipud(net) nan(nr,1); nan(1,nc+1)]);
    set(h,'edgecolor', 'none')
    set(gca, 'fontsize', 15)
    axis('off')
    colorbar

endfunction
