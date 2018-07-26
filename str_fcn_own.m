%sys x; 
D_infinite = @(x) 6.88*(x/r0).^(5/3); 
fplot(D_infinite);
xlim([0 1]);
ylim([0 350]);