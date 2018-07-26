% This procedure generates the FFT-based turbulence phase screen.Written by Xiang jingsong,China
% Reference 1:Jingsong Xiang,"Accurate compensation of the low-frequency components for the  FFT-based turbulent phase screen",Optics Express,2012.1,Vol.20,No.1,pp:681-687
% And the high frequency compensation have been added which are not described in Ref.1

function  Xiang_Jingsong_FFT_COV_low_high_com()
Dvar_theory_L0_fun=inline('power(L0/r0,5/3)*power(2,1/6)*gamma(11/6)*power(pi,-8/3)*power(24/5*gamma(6/5),5/6)*(gamma(5/6)*power(2,-1/6)-power(2*pi*r/L0,5/6).*besselk(5/6,2*pi*r/L0))','r','r0','L0');
Dvar_theory_inf=inline('6.8838771823*power(r/r0,5/3)','r','r0');

aspect_ratio=2;
My=128;Mx=aspect_ratio*My;
LMx=aspect_ratio*8+1;LMy=1*8+1;

if Mx==My
    N_set_zero_x=3;
else 
    N_set_zero_x=3*Mx/My+1;
end
N_set_zero_y=3;

Dy=1;
Dx=Dy*Mx/My;
dxy=Dx/Mx;

r0=0.2;
L0=100;
l0=0.0001;

fprintf('Dx=%6.2f,Dy=%6.2f,L0=%6.2f,r0=%6.3f,l0=%6.3e,dxy=%6.3e\n',Dx,Dy,L0,r0,l0,dxy); 



Mx0=2*Mx;My0=2*My;
HMx=Mx/2;HMy=My/2;
HMx0=Mx0/2;HMy0=My0/2;
HLMx=(LMx-1)/2;HLMy=(LMy-1)/2;
HN_set_zero_x=floor(N_set_zero_x/2);
HN_set_zero_y=floor(N_set_zero_y/2);
dNx=HMx/HLMx;dNy=HMy/HLMy;

Dx0=2*Dx;Dy0=2*Dy;
dkx0=2*pi/Dx0;dky0=2*pi/Dy0;


k0=2*pi/L0;
k02=k0*k0;
km=2*pi/l0;
km2=km*km;


%high frequency compensation
%Reference 3:???,"?????????????????????”??????2014.10,Vol.34?No.10,1001003
%because the power spectrum larger than pi/dxy is not contained in the stardard FFT phase screen,
%this high frequency compensation can contain the energy which frequency larger than pi/dxy 
n_high_com=6;
[txx,tyy] = meshgrid(-n_high_com:n_high_com);
rh=dxy*sqrt(txx.^2+tyy.^2);
%Integral the autocorrelation of the high frequency(larger than pi/dxy?power spectrum,
%and then through FFT to obtain the need compensated power spectrum in the lower frequency Mx0*My0 grid.  
cov_high=covfun_up_array(rh,n_high_com,n_high_com,r0,L0,l0,pi/dxy);
temp_covMM0=zeros(Mx0,My0);
temp_covMM0(HMx0+1-n_high_com:HMx0+1+n_high_com,HMy0+1-n_high_com:HMy0+1+n_high_com)=cov_high;
psd2_highcom=ifftshift(real(ifft2(ifftshift(temp_covMM0))));


% generating the interpolation coordinates
sdx=-HMx*dxy;sdy=-HMy*dxy;
edx=HMx*dxy; edy=HMy*dxy;
Ldx=dxy*dNx;Ldy=dxy*dNy;
[interpLMx,interpLMy] = meshgrid(sdy:Ldy:edy,sdx:Ldx:edx);
[interpMx,interpMy] = meshgrid(dxy*(-HMy:1:My-HMy-1),dxy*(-HMx:1:Mx-HMx-1));
 

 [xx,yy]=meshgrid(dky0*(-HMy0:1:HMy0-1),dkx0*(-HMx0:1:HMx0-1));
 mr2=xx.^2+yy.^2;
 mask_psd2=exp(-mr2/(pi/dxy)^2);
 mask_psd2(mask_psd2>=exp(-1))=1;
 mask_psd2(mask_psd2<exp(-1))=0;     
 %Introduce "mask_psd2" is to facilitate the calculating the high frequency autocorrelation cov_high0
 
% discreting turbulence power spectrum for FFT-based phase screen 
psd2=0.489837*dkx0*dky0*power(r0,-5/3)*power(mr2+k02,-11/6).*exp(-mr2/km2);
psd2=mask_psd2.*psd2+psd2_highcom;
%psd2 is power spectrum with high frequency compensation
%the central 3*3 grid of psd2 is set to zero
psd2(HMx0+1-HN_set_zero_x:HMx0-HN_set_zero_x+N_set_zero_x,HMy0+1-HN_set_zero_y:HMy0-HN_set_zero_y+N_set_zero_y)=0;

shift_psd2=ifftshift(psd2);
shift_psd=sqrt(shift_psd2);

r_Dvar_theory=dxy*[0,1,sqrt(2),2,sqrt(5),2*sqrt(2),3,4,5,6:4:1.6*sqrt(Mx^2+My^2)];
[Dphase_theory,cov00]=Dphasefun(r_Dvar_theory,r0,L0,l0);

% computing the autocorrelation of the FFT phase screen 
covMM0=Mx0*My0*real(ifft2(shift_psd2));
cov_highLMLM=covMM0(1:dNx:Mx+1,1:dNy:My+1);

[XX,YY] = meshgrid((0:LMy-1)*Ldy,(0:LMx-1)*Ldx);
ri=sqrt(XX.^2+YY.^2);
Dvar_theoryLMLM=interp1(r_Dvar_theory,Dphase_theory,ri,'spline');
cov_theoryLMLM=cov00-0.5*Dvar_theoryLMLM;
cov_lowLMLM=cov_theoryLMLM-cov_highLMLM;

% computing the generating coefficients of the low resolution low frequency phase screen 
coeff_phase_create=phase_screen_low_coeff(LMx,LMy,cov_lowLMLM,cov_theoryLMLM);

%The theoretical structure functions with finite outer scale L0 and nonzero inner scale l0 (Condition:L0>>l0)
rx=(1:Mx-1)*dxy;
ry=(1:My-1)*dxy;
Dvar_theory_x=interp1(r_Dvar_theory,Dphase_theory,rx,'spline');
Dvar_theory_y=interp1(r_Dvar_theory,Dphase_theory,ry,'spline');

fprintf('preparatory work is completed\n\n');
 
snum1=1000;
snum2=1000;
allDvarx1=0;allDvary1=0;

 for ipack=1:snum2
    time0=clock;
   [tvarx1,tvary1]=simulation_one_group(snum1,Mx0,My0,LMx,LMy,shift_psd,coeff_phase_create,interpLMx,interpLMy,interpMx,interpMy);
   
    num_sim=ipack*snum1;
    
    allDvarx1=allDvarx1+tvarx1; 
    allDvary1=allDvary1+tvary1; 
  
    Dvarx1=allDvarx1/num_sim;
    Dvary1=allDvary1/num_sim;
    
    N_x1=length(tvarx1);
    N_y1=length(tvary1);
       
    errx1=(Dvar_theory_x(1:N_x1)-Dvarx1(1:N_x1))./Dvar_theory_x(1:N_x1);
    erry1=(Dvar_theory_y(1:N_y1)-Dvary1(1:N_y1))./Dvar_theory_y(1:N_y1);
 
    maxerrx1=max(abs(errx1(1:end)));
    maxerry1=max(abs(erry1(1:end)));
  
    
    Time_group=etime(clock,time0);
    fprintf('Prog=%7d/%d, t=%6.2f, max_errx1=%6.4f, max_erry1=%6.4f\n',num_sim,snum1*snum2,Time_group,maxerrx1,maxerry1); 
 end
  
plot(rx(1:length(errx1)),errx1,'-b',ry(1:length(erry1)),erry1,'-.b');





function [coeff_phase_create,phasecov_LM_LM]=phase_screen_low_coeff(LMx,LMy,cov_lowLMLM,cov_theoryLMLM)
%generating autocorrelation of the low resolution phase screen
time0=clock;
Dvar_lowLMLM=2*(cov_lowLMLM(1,1)-cov_lowLMLM);
Dvar_theoryLMLM=2*(cov_theoryLMLM(1,1)-cov_theoryLMLM);

phasecov_LM_LM=zeros(LMx*LMy,LMx*LMy);
phasecov_all_LM_LM=zeros(LMx*LMy,LMx*LMy);
for ix=1:LMx
    for iy=1:LMy
      imm=(iy-1)*LMx+ix;
      for iix=1:LMx
        for iiy=1:LMy
         iimm=(iiy-1)*LMx+iix; 
          phasecov_LM_LM(imm,iimm)=cov_lowLMLM(1+abs(ix-iix),1+abs(iy-iiy));
          phasecov_all_LM_LM(imm,iimm)=cov_theoryLMLM(1+abs(ix-iix),1+abs(iy-iiy));
%          phasecov_LM_LM(imm,iimm)=0.5*(Dvar_lowLMLM(ix,iy)+Dvar_lowLMLM(iix,iiy)-Dvar_lowLMLM(1+abs(ix-iix),1+abs(iy-iiy)));
%          phasecov_all_LM_LM(imm,iimm)=0.5*(Dvar_theoryLMLM(ix,iy)+Dvar_theoryLMLM(iix,iiy)-Dvar_theoryLMLM(1+abs(ix-iix),1+abs(iy-iiy)));
        end
     end
 end
end
[V,D]=eig(phasecov_LM_LM);
coeff_phase_create=(V)*real(sqrt(D));


Dvar_creat=zeros(LMx*LMy,LMx*LMy);
Dvar_original=zeros(LMx*LMy,LMx*LMy);
Dvar_theory=zeros(LMx*LMy,LMx*LMy);
phasecov_creat=(coeff_phase_create)*transpose(coeff_phase_create);
for ix=1:LMx*LMy
    for iy=1:LMx*LMy  
     Dvar_creat(ix,iy)=phasecov_creat(ix,ix)+phasecov_creat(iy,iy)-2*phasecov_creat(ix,iy);
     Dvar_original(ix,iy)=phasecov_LM_LM(ix,ix)+phasecov_LM_LM(iy,iy)-2*phasecov_LM_LM(ix,iy);
     Dvar_theory(ix,iy)=phasecov_all_LM_LM(ix,ix)+phasecov_all_LM_LM(iy,iy)-2*phasecov_all_LM_LM(ix,iy);
    end
end
terr=abs((Dvar_original-Dvar_creat)./(Dvar_theory+1.0e-20));
terr=terr.*(1-eye(LMx*LMy));
errSVD=max(max(terr));
fprintf('Time of SVD =%6.2f, Minimum eigenvalue =%6.4f, Maxerr =%6.4e\n',etime(clock,time0),min(diag(D)),errSVD); 




function [Dvarx1,Dvary1]=simulation_one_group(num,Mx0,My0,LMx,LMy,shift_psd,coeff_phase_create,interpLMx,interpLMy,interpMx,interpMy)
 
  Mx=Mx0/2;My=My0/2;
  HMx=Mx/2;HMy=My/2;
  dNx=Mx/(LMx-1);dNy=My/(LMy-1);
  
  startx1=2;starty1=2;
  
  Dvarx1=zeros(1,Mx-startx1);
  Dvary1=zeros(1,My-starty1);
  
  xjs1=clock;
  for inum=1:2:num
           
    randomfft=randn(Mx0,My0)+i*randn(Mx0,My0);
    Pt=fft2(shift_psd.*randomfft);
    P_FFT=Pt(HMx+(1:Mx),HMy+(1:My)); 

    %generating the low frequency compensating phase screen based on covariance method
    randomcov=randn(LMx*LMy,1)+i*randn(LMx*LMy,1);
    LP_one=(coeff_phase_create*randomcov);
%     transform column matrix of (1,LM*LM) to square matrix of (LM,LM) 
    LP=reshape(LP_one,LMx,LMy);
   %generating high resolution compensating phase screen using interpolation method
    P_add=interp2(interpLMx,interpLMy,LP,interpMx,interpMy,'spline'); 
%     the final phase screen
    P=P_FFT+P_add;
    
    P1=real(P);
    P2=imag(P);
% 
    tx1=(P1(startx1,starty1)-P1(startx1+1:Mx,starty1)).^2;
    tx2=(P2(startx1,starty1)-P2(startx1+1:Mx,starty1)).^2;
    Dvarx1=Dvarx1+transpose(tx1+tx2);
    
    ty1=(P1(startx1,starty1)-P1(startx1,starty1+1:My)).^2;
    ty2=(P2(startx1,starty1)-P2(startx1,starty1+1:My)).^2;
    Dvary1=Dvary1+ty1+ty2;

  end

  
  
function cov_theory=covfun_up_array(druo,HNx,HNy,r0,L0,l0,kdown)
covfun=inline('besselj(0,r*t).*power(t.^2+k02,-11/6).*exp(-t.^2/km2).*t','t','r','k02','km2');
 time0=clock;

 HN=HNx;
 N=2*HN+1;
 
 k0=2*pi/L0;
 km=2*pi/l0;
 k02=k0^2;
 km2=km^2;
 xs=6.1554729787*power(r0,-5/3);
   
 scale=5; num=2000;
 num_iter=floor(log10(6*km/kdown)/log10(scale));
 
 rx=druo(HN+1,HN+1:N);rx(1)=1.0e-20;
 
 cov_theory=zeros(N,N);
for ix=HN+1:N
    for iy=HN+1:ix
    r=druo(ix,iy);
    st=kdown;et=kdown*scale;dx=(et-st)/num;
    x=(st+dx/2.0:dx:et);
    tall=sum(covfun(x,r,k02,km2))*dx;
    for in=1:num_iter
       st=et; et=et*scale; 
       dx=(et-st)/num;
       x=(st+dx/2.0:dx:et);
       tall=tall+sum(covfun(x,r,k02,km2).*dx);   
    end 
    cov_theory(ix,iy)=0.5*xs*tall;
    end
end

for ix=HN+1:N
    for iy=ix:N
     cov_theory(ix,iy)=cov_theory(iy,ix);
    end
end
for ix=1:HN
    for iy=HN+1:N
     cov_theory(ix,iy)=cov_theory(N+1-ix,iy);
    end
end
for ix=1:N
    for iy=1:HN
     cov_theory(ix,iy)=cov_theory(ix,N+1-iy);
    end
end

function [Dphase_theory,cov00]=Dphasefun(druo,r0,L0,l0)
 Dpfun=inline('(1-besselj(0,r*t)).*power(t.^2+k02,-11/6).*exp(-t.^2/km2).*t','t','r','k02','km2');
 Dp00=inline('power(t.^2+k02,-11/6).*exp(-t.^2/km2).*t','t','k02','km2');
 xjs0=clock;

 k0=2*pi/L0;
 km=2*pi/l0;
 k02=k0^2;
 km2=km^2;
 % 2*power(24/5*gamma(6/5),5/6)=6.8838771823
 % 0.033*8*pi^3*power(2,-8/3)*power(24/5*gamma(6/5),-5/6)/gamma(11/6)^2=0.42329413
 % 0.033*2*pi/0.42329413=0.489837;
 % 0.033*8*pi*pi/0.42329413=6.155473;
 xs=6.155473*power(r0,-5/3);
    
 scale=5;
 num=5000;
 num_iter=1+floor(log10(6*km/k0)/log10(scale));
 
 st=0;et=k0;
 dt=(et-st)/num;
 x=(st+dt/2:dt:et);
 tall=sum(Dp00(x,k02,km2))*dt;
 for i=1:num_iter
       st=et; et=et*scale;  
       dt=(et-st)/num;
       x=(st+dt/2:dt:et);
       t=sum(Dp00(x,k02,km2))*dt;
       tall=tall+t;
 end 
  cov00=0.5*xs*tall;

 
for ip=1:length(druo) 
   r=druo(ip);
 if r>1.0e-12
   st=0;et=k0;
   dt=(et-st)/num;
   x=(st+dt/2:dt:et);
   tall=sum(Dpfun(x,r,k02,km2))*dt;
   for i=1:num_iter
       st=et; et=et*scale;  
       dt=(et-st)/num;
       x=(st+dt/2:dt:et);
       t=sum(Dpfun(x,r,k02,km2))*dt;
       tall=tall+t;
   end 
    Dphase_theory(ip)=xs*tall;
 else
    Dphase_theory(ip)=0;  
 end
end
% cov_theory=cov00-0.5*Dphase_theory;
 