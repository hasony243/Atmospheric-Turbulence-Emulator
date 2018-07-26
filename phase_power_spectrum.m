%plot phase power spectrum to proof the limatation of FFT method 

p3 = fplot(@(x) 0.033*((x^2+L0^(-2))^(-11/6)), [0 1.5], 'm');
legend('L0=100')
ylim([0 500]); 
hold on; 
grid on;
xlabel('wavenumber');
ylabel('phase power spectrum'); 
