D = 2; % length of one side of square phase screen [m]
r0 = 0.1; % coherence diameter [m]
N = 256; % number of grid points per side
L0 = 100; % outer scale [m]
l0 = 0.01;% inner scale [m]

delta = D/N; % grid spacing [m]

% spatial grid
x = (-N/2 : N/2-1) * delta;
y = x;

%[x y] = meshgrid(x);
%w=2; % width of rect
%mask=ones(N);
%A=rect(x/w).*rect(y/w); %signal

% generate a random draw of an atmospheric phase screen
B_corr_ft = twoD_autocor_ft(r0, N, delta, L0, l0);
