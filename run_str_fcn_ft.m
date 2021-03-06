L = 2; % length of one side of square phase screen [m]
r0 = 0.1; % coherence diameter [m]
N = 256; % number of grid points per side
L0 = 100; % outer scale [m]
l0 = 0.1;% inner scale [m]
%L is grid size
delta = L/N; % grid spacing [m]
del_f = 1/(N*delta);
% spatial gridsuqat
x = (-N/2 : N/2-1) * delta;
y = x;

%[x y] = meshgrid(x);
%w=2; % width of rect
%mask=ones(N);
%A=rect(x/w).*rect(y/w); %signal

% generate a random draw of an atmospheric phase screen
%for i = 1:30
B_corr_ft = corr_ft(r0, N, delta, L0, l0, L);
B_corr_ft_0 = 1; 
D = 2*(B_corr_ft_0-B_corr_ft); 

%figure 
%mesh(D); 
%t = x - y;
%mesh(x,y,D,t); 
%colormap(line(6)); 
%tol = 1e-8; 
%mask=abs(t)<tol; 
%x2 = x(mask); 
%y2 = y(mask); 
%z2 = D(mask); 
%plot3(x2,y2,z2); 
%plot(x2,z2/r0); 

%D_min=D-min(D);

str = D(N/2+1,N/2+1:N);
x = (1:N/2)*delta;
str = str+abs(min(str));
plot(x,str); xlim([0 1]);
%str=D(N/2+1,N/2+1:N);
%x = (N/2+1:N)*delta; 
%plot(x,str);
%M=min(str);
%str=str-M; 
%plot(x,str);
%xlim([0 1]);
%w2 = ift2(W.*conj(W), delta_f);
%mesh(B_corr_ft) 
%t = x-y;
%figure
%mesh(x,y,B_corr_ft,t);
%colormap(lines(6));
%tol=1e-8; 
%mask=abs(t)<tol; 
%x2=x(mask);
%y2=y(mask);
%z2=C(mask);
%figure
%plot3(x2,y2,z2)
%figure
%plot(x2,z2)
%xlim([0 2]);
%hold on;
%p3 = fplot(@(x) 6.88*((x/r0)^(5/3)), [0 1], 'r', 'DisplayName','Theory');
%legend;
%ylim([0 350]);