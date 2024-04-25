%% CRLB as a function of the SNR (for different NR values) for the three methods
clear all
addpath('/Users/carla.macias/Desktop/ris-crlb-optimization/src/optimization_functions');
warning('off','all')
xBS=0; yBS=0; xRIS=5; yRIS=5; xM_real=10; yM_real=0.5;
Nr=51; Nt= 21;

f = 28; j=sqrt(-1); Q=0.2;
c = physconst('LightSpeed'); lambda = c/(f*10^9);
Lbr=1; Lrm=1; delta=lambda/2; mu=3; B=(Nr-1)/2; M=(Nt-1)/2;
B_vect=-B:B; M_vect=-M:M;

% Distance between BS and scatters & between scatteres and RIS
r_br=sqrt((xRIS-xBS)^2+(yRIS-yBS)^2); d_br=r_br;
% Angle of departure between RIS and BS
AODbr_grad=acosd((xRIS-xBS)/r_br)+90; AODbr=AODbr_grad*(2*pi)/360;
% Angle of arrival between RIS and BS (deg)
AOAbr_grad=180+AODbr_grad;  AOAbr=AOAbr_grad*(2*pi)/360;

% Distance between RIS and scatters & scatteres an UE
r_rm=sqrt((xM_real-xRIS)^2+(yM_real-yRIS)^2); d_rm= r_rm;
% Angle of departure between UE and RIS (deg)
AODrm_grad=acosd((abs(yRIS-yM_real))/r_rm);  AODrm=AODrm_grad*(2*pi)/360;
% Angle of arrival between UE and RIS (deg)
AOArm_grad=180+AODrm_grad;  AOArm=AOArm_grad*(2*pi)/360;


% Cota en funcio de la SNR -grafica-
%Cota CRLB
  CRLB_vect_tot=[];
for Nris=[10, 30, 50, 80, 100, 150, 200]
    CRLB_vect=[];
   
    w=1;fb_vect=nan(1,length(B_vect));
    for bsum1=-B:B
        fb_vect(w)=bsum1*omega(AODbr,delta,lambda)+(bsum1^2)*gamma1(AODbr,r_br,delta,lambda);
        w=w+1;
    end
    Fb_sum=sum(fb_vect);

    w=1;fm_vect=nan(1,length(M_vect));
    for msum1=-M:M
        fm_vect(w)=msum1*omega(AOArm,delta,lambda)+(msum1^2)*gamma1(AOArm,d_rm,delta,lambda);
        w=w+1;
    end
    Tm_sum=sum(fm_vect);

    phasesRIS_SNR=nan(1,Nris);
    for r=1:Nris
        Gr=r*omega(AOAbr,delta,lambda)+r^(2)*gamma1(AOAbr,r_br,delta,lambda)+...
            r*omega(AODrm,delta,lambda)+r^(2)*gamma1(AODrm,d_rm,delta,lambda);
        phasesRIS_SNR(r)=((2*M+1)*Fb_sum + (2*B+1)*(2*M+1)*Gr + (2*B+1)*Tm_sum)*(((2*M+1)*(2*B+1))^(-1));
    end


    for SNR=-20:5:20

        func=@(phaseRIS)CRLBfunction(xBS,yBS,xRIS,yRIS,xM_real,yM_real,phaseRIS,Nris,SNR, Q);
        [fases,fval,exitflag,output]=fminsearch(func,phasesRIS_SNR);
        CRLB_vect=[CRLB_vect sqrt(fval)];
    end
    CRLB_vect_tot=[CRLB_vect_tot; CRLB_vect];
end
%%
fig=figure(1);
semilogy(-20:5:20,CRLB_vect_tot(1,:),'-o')
hold on
semilogy(-20:5:20,CRLB_vect_tot(2,:),'-o')
hold on
semilogy(-20:5:20,CRLB_vect_tot(3,:),'-o')
hold on
semilogy(-20:5:20,CRLB_vect_tot(4,:),'-o')
hold on
semilogy(-20:5:20,CRLB_vect_tot(5,:),'-o')
hold on
grid on
legend('Nris=10','Nris=30','Nris=50','Nris=80','Nris=100')
ylabel('CRB')
xlabel('SNR [dB]')


save('SNRvsCRLB')
