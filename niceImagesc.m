function niceImagesc (net, bluemap)

    figure;
    [nr nc] = size(net);
    net(net == 0) = nan;    %pcolor can be set to ignore nans
    b = linspace(0.1, 0.1 , 50);
    bb = linspace(0.1, 0.8, 25);
    b = [b bb];
    b = flipud(b);
    bmap = [b; b; [b(1:30) ones(1, 45)]];
    
    if(bluemap)
        colormap(bmap');
    end
    
    h = pcolor([flipud(net) nan(nr,1); nan(1,nc+1)]);
    set(h,'edgecolor', 'none')
    set(gca, 'fontsize', 15)
    axis('off')
%    colorbar

endfunction
