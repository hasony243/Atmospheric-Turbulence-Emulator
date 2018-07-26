load verg1; 
Ts = 0.005; %sample interval = 5ms 
[nu, N] = size(data_out);

t = [1:N]*Ts; %generate time vector 

plot(t, data_out, 'k'); 
hold on; 

avg = mean(data_out);
avg(1:10);

plot(t, avg-3, 'k'); 
xlabel('Time(sec)'); 
ylabel('Eye position');

plot([.43 .43], [0 5], 'k');
text(1, 1.2, 'Average Data');