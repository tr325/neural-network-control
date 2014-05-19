%## Copyright (C) 2014 Tom Renner

% Takes input of stabilised W matrix (from c++ code), and simulates network dynamics
function plotDynamicEvolution(W)
    
    W = squeeze(W);
    [dim p] = size(W);
    A = W - eye(dim);
    
    %to get largest direction of amplification (see G's paper)
    Q = lyap(A', 2*eye(dim));
    [a e] = eigs(Q, 1, 'lm');

%    a = normrnd(0, 0.5, dim, 1);        %input preparatory signal
    x = zeros(dim, 1);
    xIP = a - W*a;
    
    X = [];         %container for time-evolving x vectors
    t = 0.01;       % time step
    

    %implement eqn xdot = W*x; 
    %while input is held
    for i = 1:1000
        if i >500
            X = [X x];
        end
        xdot = A*x + xIP;
        x = x + t*xdot;
    end
    %after input is released
    for i = 1:999
        X = [X x];
        xdot = A*x;     
        x = x + t*xdot;
    end
        
    
    figure
    hold all
    % only 10 for clarity of image
    for i = 25:35
        plot(X(i, :))
    end
    axis([0 1500])
    xlabel('Time (arbitrary units)', 'fontsize', 15);
%    ylabel('Value of elements of x vector', 'fontsize', 15);
    hold off

%    figure
%    plot(norm(X,2,'cols').*sign(sum(X)), '.')
%    
    
end