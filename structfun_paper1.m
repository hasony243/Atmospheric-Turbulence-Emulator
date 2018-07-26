N_1 = 256;
for F=1:10:511
    d1=0;
    n=0;
    for a=1:N_1-F
        for b=1:N_1-F
            n=n+2;
            d1=d1+(phz(a+F,b)-phz(a,b))^2;
            d1=d1+(phz(a,b+F)-phz(a,b))^2;
        end
    end
    D1(1,(F+9)/10)=d1/n;
end

D = D1; 
r = (1:length(D)); 
plot(r-1,D)