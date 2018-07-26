%example_strfcn2_rect.m 

N = 256; %number of samples  
L = 16;  %grid size [m]
delta = L/N; %sample spacing [m]
F = 1/L; %frequency-domain grid spacing[1/m]
x = [-N/2 : N/2-1]*delta; 
[x y] = meshgrid(x);
w = 2; %width of rectangle
%A = rect(x/2).*rect(y/w);
A = lambdaWrapped;
%A = phz;
mask = ones(N); 
%perform digital structure function 
C = str_fcn2_ft(A, mask, delta);
C = real(C);
%C = str_fcn2_ft(A, mask, delta)/power(delta,2);
%continuous structure function 
%C_cont = 2 * w^2 * (1 - tri(x/w) .* tri(y/w));
figure
mesh(C) 
t = x-y;
figure
mesh(x,y,C,t);
colormap(lines(6));
tol=1e-8; 
mask=abs(t)<tol; 
x2=x(mask);
y2=y(mask);
z2=C(mask);
figure
plot3(x2,y2,z2)
figure
plot(x2,z2/r0)
xlim([0 1]);
hold on;
p3 = fplot(@(x) 6.88*((x/r0)^(5/3)), [0 1], 'r', 'DisplayName','Theory');
