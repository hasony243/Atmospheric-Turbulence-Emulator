%clear all; 
D_fun = zeros(1,round(sqrt(2)*N));
for i = 1:1000
    %D_fun = zeros(1,256);
    run_ft_sh;
    %run_ft;
    structfun_def; %got D from this line
    %structfun_def_v1; %got D from this line
    D_fun = D_fun + D; 
end 
D = D_fun ./ 1000;  
r = (1:length(D_fun)) * delta;
figure; 
plot(r,D); xlim([0 1]);
%plot(D)