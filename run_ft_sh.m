% example_ft_phase_screen.m

D = 2; % length of one side of square phase screen [m]
r0 = 0.1; % coherence diameter [m]
N = 256; % number of grid points per side
L0 = 100; % outer scale [m]
l0 = 0.01;% inner scale [m]

delta = D/N; % grid spacing [m]
% spatial grid
x = (-N/2 : N/2-1) * delta;
y = x;
% generate a random draw of an atmospheric phase screen
[phz_lo phz_hi] = ft_sh_phase_screen(r0, N, delta, L0, l0);
phz_sub = phz_lo + phz_hi;

%map = parula(128);

%angleInradians=degtorad(phz_sub);

%colormap(jet(64));
lambdaWrapped = wrapToPi(phz_sub);
%x = [-D/2 D/2]; 
%y = [-D/2 D/2];
%clims=[0 2];
%imagesc(phz_sub)
imagesc(lambdaWrapped)
%imagesc(angleInradians);
%figure