%%plot the structure function for the case of a finite outer scale 
%r0 = 0.1; 
t1 = 6.16*(r0^(-5/3)); 
%t1 = power(L0/r0,5/3);
t2 = (3/5)*((L0/(2*pi))^(5/3));
%t2 = power(2,1/6)*gamma(11/6)*power(pi,-8/3); 
t3 = ((L0/(4*pi))^(5/6))/gamma(11/6); 
%t3 = power(24/5*gamma(6/5),5/6);
%t4 = gamma(5/6)*power(2,-1/6);
%t5 = power(2*pi*x/L0,5/6).*besselk(5/6, 2*pi*x/L0);
%syms x;

p1 = fplot(@(x) t1*(t2-(x^(5/6))*t3*besselk((5/6),(2*pi*x/L0))), 'm', [0 1]);  
%p1 = fplot(@(x) t1*t2*t3*(t4-t5));

legend('D-theory proposed by [1] with L0 = 100m')
%fplot(@(r) power(L0/r0,5/3)*power(2,1/6)*gamma(11/6)*power(pi,-8/3)*power(24/5*gamma(6/5),5/6)*...
%(gamma(5/6)*power(2,-1/6)-power(2*pi*r/L0,5/6).*besselk(5/6,2*pi*r/L0)));
%hold on;
%p3 = fplot(@(x) 6.88*(x^(5/3)), [0 13], 'g');
xlabel('r');
ylabel('D(r)');
%ylim([0 350]);


%%plot the structure function for the case of a finite outer scale 
%t1 = 6.16*(r0^(5/3)); 
%t2 = (3/5)*((L0/(2*pi))^(5/3));
%t3 = ((L0*r0/(4*pi))^(5/6))/gamma(11/6); 
%syms x;

%p2 = fplot(@(x) t1*(t2-(x^(5/6))*t3*besselk((5/6),(2*pi*x/L0))), 'r', [0 13]);  
%ylim([0 400]);
%xlabel('r');
%ylabel('D(r)');

%Dvar_theory_L0_fun=inline('power(L0/r0,5/3)*power(2,1/6
%)*gamma(11/6)*power(pi,-8/3)*power(24/5*gamma(6/5),5/6)*(
%gamma(5/6)*power(2,-1/6)-power(2*pi*r/L0,5/6).*besselk(5/6,2*pi*r/L0))','r','r0','L0');

%fplot(@(x) power(L0/r0,5/3)*power(2,1/6)*gamma(11/6)/power(pi,11/6)*power(24/5*gamma(6/5),5/6))...
 %   *(gamma(6/5)/power(2,1/6)-power((2*pix/L0),5/6)*besselk(5/6,2*pi*x/L0));
 
