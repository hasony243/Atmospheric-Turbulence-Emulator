phi = lambdaWrapped; 
%phi = phz;
N = size(phi,1); % number of samples  
L = 2;  % grid size [m]
delta = L/N; % sample spacing [m]
coords = randi(N,1,2); %generate the coordinates of specific point 
friedPar = 1000213012; %min_dist
a = 1; b = 1; 
for i = 1:N
   for j = 1:N
    if abs(phi(coords(1,1),coords(1,2)) - phi(i,j)) <= 1 && ...
            abs(phi(coords(1,1),coords(1,2)) - phi(i,j)) >= 0.95 && i~=a && j~=b 
        if sqrt((coords(1,1)-i)^2+(coords(1,2)-j)^2)*delta >= friedPar
            continue
        else
            a = i; b = j;
            friedPar = sqrt((coords(1,1)-a)^2+(coords(1,2)-b)^2)*delta;
        end
     end
    end   
end
friedPar;