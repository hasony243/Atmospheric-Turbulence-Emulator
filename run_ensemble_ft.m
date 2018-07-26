%run_ensemble_fourier 
tic

%D = 500; % length of one side of square phase screen [m]
%r0 = 0.1; % coherence diameter [m]
%N = 256; % number of grid points per side
%L0 = 100; % outer scale [m]
%l0 = 0.01;% inner scale [m]

%delta = D/N; % grid spacing [m]

% spatial gridsuqat
%x = (-N/2 : N/2-1) * delta;
%y = -x;
i = 1;
for i = 1:i
    clear all; 
    D = 2; % length of one side of square phase screen [m]
    r0 = 0.1; % coherence diameter [m]
    N = 256; % number of grid points per side
    L0 = 100; % outer scale [m]
    l0 = 0.01;% inner scale [m]

    delta = D/N; % grid spacing [m]

% spatial gridsuqat
    x = (-N/2 : N/2-1) * delta;
    y = -x;
    
    phz = ft_phase_screen(r0, N, delta, L0, l0); 
    
    phz_total = 0;
    
    phz_total = phz_total + phz; 
end

mask = ones(N);
delta=D/N;
phz_ensemble=phz_total/i;

%display the ensemble phase screen
colormap(jet(64));
lambdaWrapped = wrapToPi(phz_ensemble);
%for j = 1:N
%    for i = 1:N
%        lambdaWrapped_prime(i,j)=lambdaWrapped(N+1-i,j);
%    end 
%end
clims=[-pi pi];
%phz_image=imagesc(lambdaWrapped,clims);
%C = str_fcn2_ft(phz_ensemble, mask, delta);

%D = real(C);
%D_sub = real(C_sub);

%p1=plot(x/(r0), D, 'g--');
%xlim([0 13])
%ylim([0 10])
toc
