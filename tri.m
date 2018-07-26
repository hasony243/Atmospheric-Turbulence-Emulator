function y = tri(t)
% function y = tri(t)
    t = abs(t);
    y = zeros(size(t));
    idx = find(t < 1.0); 
    y(idx) = 1.0 - t(idx);