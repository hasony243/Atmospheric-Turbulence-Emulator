%calculation phase structure function by definition 
phi = phz_sub; 
%phi = lambdaWrapped;
N = size(phi,1); % number of samples  
L = 2;  % grid size [m]
delta = L/N; % sample spacing [m]
%D = zeros(1,N);
D = zeros(1,round(sqrt(2)*N)); % output function
count = D; % for computing mean

for n = 1:1e6 % find a good amount of points here, the more points the better the estimate
   coords = randi(N,2,2);
   r = round(norm(coords(1,:) - coords(2,:)));
   if r<1
      continue % skip if the two coordinates are the same
   end
   d = phi(coords(1,1),coords(1,2)) - phi(coords(2,1),coords(2,2));
   %d = mod(abs(d),pi); % might not need this, depending on how A is constructed
   D(r) = D(r) + d.^2;
   %D(r) = D(r).^2; 
   count(r) = count(r) + 1;
end
I = count > 0;
D(I) = D(I) ./ count(I); % do not divide by 0, some bins might not have any samples
I = count < 100;
D(I) = 0; % ignore poor estimates

r = (1:length(D)) * delta; 
plot(r,D); 
xlim([0 1]);