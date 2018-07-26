%tic
D = zeros(1,256);
%ph = phz;
ph = lambdaWrapped;
for delta_r = 1:256 
    % a, b is position of x(a,b) 
    for a = 1:256 
        for b = 1:256 
            count = 0; 
            D1_delta_r = 0;
            D2_delta_r = 0;
            for i = (a-delta_r):(a+delta_r) 
                if (i <= 0)||(i > 256)
                    continue
                end
                if (a-delta_r>0) 
                    D1_delta_r = power(ph(a,b)-ph(a-delta_r,i),2) + D1_delta_r; 
                    count = count + 1;
                end 
                if (a+delta_r<256)
                    D2_delta_r = power(ph(a,b)-ph(a+delta_r,i),2) + D2_delta_r;
                    count = count + 1;
                end 
                %end 
            end
            
            D3_delta_r = 0;
            D4_delta_r = 0;
            for i = (b-delta_r+1):(b+delta_r-1)
                if ( i <= 0)||( i>256 )
                    continue
                end 
                if (b-delta_r>0)
                    D3_delta_r = power(ph(a,b)-ph(i,b-delta_r),2) + D3_delta_r;
                    count = count + 1;
                end 
                if (b+delta_r<256)
                    D4_delta_r = power(ph(a,b)-ph(i,b+delta_r),2) + D4_delta_r;
                    count = count + 1;
                end 
            end
            %end
        end 
    end 
    D(1,delta_r) = (D1_delta_r+D2_delta_r+D3_delta_r+D4_delta_r)/count;
end
figure
plot(D)
%toc 

