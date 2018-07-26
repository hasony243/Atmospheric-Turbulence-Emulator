%strfcn2_rect.m 

%N = 512; 
L = 100; 
delta = L/N; 
F = 1/L; 
A = phz;
A_prime = lambdaWrapped;
%A_sub = phz_sub;
mask = ones(N); 
x = [-N/2 : N/2-1]*delta;
C = str_fcn2_ft(A, mask, delta)/(delta^2);
%C_sub = str_fcn2_ft(A_sub, mask, delta);

D = real(C);
%D_sub = real(C_sub);
row_center = D(N/2+1,:); 
row_center_plot=row_center(1,N/2+1:end);

p1= plot(x, row_center,'- o');
xlim([0 1]);
ylim([0 500]);
%p1=plot(x/(r0), D);
%plot(x)
legend('Simulated FFT method - N=256 - D=100')

%hold on;
%p2=plot(x/r0,D_sub, 'm*'); 
%legend('Simulated SH')
hold on;

%x = [0:10];
c;
%legend('Theory')
xlabel('r');
ylabel('D_\Theta');
%p3=plot(D_theory);
hold off;
%legend('Simulated FT','Theory');
%figure;
%p4 = plot(y,x);

%xlim([0 13])
%ylim([0 500])
grid on;
%hold on;
%legend('Simulated FT','Simulated SH','Theory');


