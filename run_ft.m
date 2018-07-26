%clear all;
%tic

% example_ft_phase_screen.m

D = 2; % length of one side of square phase screen [m]
r0 = 0.1; % coherence diameter [m]
N = 256; % number of grid points per side
L0 = 100; % outer scale [m]
l0 = 0.01;% inner scale [m]

delta = D/N; % grid spacing [m]
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
phz = ft_phase_screen(r0, N, delta, L0, l0);
%end
%colormap(jet(64));
%lambdaWrapped = wrapToPi(phz);
%clims=[-1 1];
%imagesc(lambdaWrapped,clims)
%B_corr_ft = ft_phase_screen(r0, N, delta, L0, l0);
%colormap(jet(64));
lambdaWrapped = wrapToPi(phz);
for j = 1:N
    for i = 1:N
        lambdaWrapped_prime(i,j)=lambdaWrapped(N+1-i,j);
    end 
end
%clims=[-pi pi];
imagesc(lambdaWrapped);
figure; 
contour(lambdaWrapped_prime,'ShowText','on');

%toc 
%splot(C)
%for i=1:40
%    imagesc(phz);
%end 
%C=str_fcn2_ft(A, phz, delta) / delta^2;




%create a temporary working folder to store the image sequency
%conworkingDir = template; 
%mkdir(workingDir)
%mk(workingDir, 'images')

%create a VideoReader to use for reading frames from the file
%shuttleVideo = VideoReader('shuttle.avi');

%lambdaWrapped = wrapToPi(phz);
%clims = [-1 1];
%imagePhase = imagesc(lambdaWrapped,clims);
%contourf(peaks);

%i1=1; 

%while hasFrame(shuttleVideo)
    %img=readFrame(shuttleVideo);
    %filename = [sprint('%03d',i1) '.fig'];
    %fullname = fullfile(workingDir,'images',filename);
    %imwrite(img,fullname);  % write out to a fig file igm1.fig ...
    %i1 = i1+1;
%end

%imageNames = dir(fullfile(workingDir,'images','*.fig'));
%imageNames = {imageNames.name}';

%outputVideo = VideoWriter(fullfile(workingDir,'shuttle_out.avi'));
%outputVideo.FrameRate = shuttleVideo.FrameRate; 
%open(ouputVideo)

%for i1=1:length(imageNames)
    %img = imread(fullfile(workingDir,'iamges',imageNames{i1}));
    %writeVideo(outputVideo, img);
%end

%clode(outputVideo);

%shuttleAvi = VideoReader(fullfile(workingDir,'shuttle_out.avi'));

%i1=1; 
%while hasFrame(shuttleAvi)
%    mov(ii) = im2frame(readFrame(shuttleAvi));
%    i1=i1+1;
%end  
%figure 
%imshow(mov(1).cdata,'Border', 'tight'); 

%movie(mov, 1, shuttleAvi.FrameRate)