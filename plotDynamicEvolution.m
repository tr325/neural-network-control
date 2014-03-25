%## Copyright (C) 2014 Tom Renner

% Takes input of stabilised W matrix (from c++ code), and simulates network dynamics
function plotDynamicEvolution(W)
    
    W = squeeze(W);
    [dim p] = size(W);
    x = normrnd(0, 1, dim, 1);
    X = [];        %container for time-evolving x vectors
    t = 0.01;       % time step
    W = W - eye(dim);

    %implement eqn xdot = W*x; 
    for i = 1:1000
        X = [X x];
        xdot = W*x;
        x = x + t*xdot;
    end

    figure
    hold all
    for i = 1:dim
        plot(X(i, :))
    end
    xlabel('Time (iterations of dynamics loop)');
    ylabel('Value of elements of x vector');
%    print -deps dyn50example.eps
    hold off

%    figure
%    plot(norm(X,2,'cols').*sign(sum(X)), '.')
%    
    
end