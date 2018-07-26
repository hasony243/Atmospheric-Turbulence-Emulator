phi = lambdaWrapped; 
%phi = phz;
N = size(phi,1); % number of samples  
L = 2;  % grid size [m]
delta = L/N; % sample spacing [m]
coords = randi(N,1,2); %generate the coordinates of specific point 
friedPar = 0; %min_dist
count = 0; 
tmpx = 0; tmpy = 0;
%x = 0; y = 0; 

for i = -N:N
    for j = -N:N 
        if (i == 0 && j ==0) 
            continue; 
        end 
        if (gcd(abs(i),abs(j))~=1) 
            continue; 
        end
        tmpx = coords(1,1) + i; tmpy = coords(1,2) + j; 
        savex = 0; savey = 0; 
        while (1 <= tmpx && tmpx <= N && 1 <= tmpy && tmpy <= N)
            if (abs(phi(tmpx,tmpy)-phi(coords(1,1),coords(1,2))) <= 1.001 && ...
                    abs(phi(tmpx, tmpy)-phi(coords(1,1),coords(1,2))) >= 0.999)
                savex = tmpx; 
                savey = tmpy; 
                break;
            end 
            tmpx = tmpx + i; 
            tmpy = tmpy + j;
        end 
        if (savex ~= 0 && savey ~= 0)
            friedPar = sqrt((savex-coords(1,1))*(savex-coords(1,1))+(savey-coords(1,2))*(savey-coords(1,2)))*delta+friedPar;
            count = count + 1; 
        end 
    end 
end 

friedPar = friedPar/count; 