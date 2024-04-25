
% For specific parameters this function gives you the CRLB
function CRLB=CRLBfunction(xBS,yBS,xRIS,yRIS,xM,yM,phaseRIS,Nris, SNR, Q, M0)
% F=sqrt(Q/2)*(randn(1)+j*randn(1));
f = 28; c = physconst('LightSpeed'); j = sqrt(-1); lambda = c/(f*10^9);
Lbr=1; Lrm=1; Nr=51;  Nt= 21; delta=lambda/2;
B=(Nr-1)/2; M=(Nt-1)/2; R=(Nris-1)/2; mu=3; 
F=1;
%% REAL POSITION

% Distance between BS and scatters & between scatteres and RIS
r_br=sqrt((xRIS-xBS)^2+(yRIS-yBS)^2); d_br=r_br;
% Angle of departure between RIS and BS
AODbr_grad=acosd((xRIS-xBS)/r_br)+90; AODbr=AODbr_grad*(2*pi)/360;
% Angle of arrival between RIS and BS (deg)
AOAbr_grad=180+AODbr_grad;  AOAbr=AOAbr_grad*(2*pi)/360;

% Distance between RIS and scatters & scatteres an UE
r_rm=sqrt((xM-xRIS)^2+(yM-yRIS)^2); d_rm= r_rm;
% Angle of departure between UE and RIS (deg)
AODrm_grad=acosd((abs(yRIS-yM))/r_rm);  AODrm=AODrm_grad*(2*pi)/360;
% Angle of arrival between UE and RIS (deg)
AOArm_grad=180+AODrm_grad;  AOArm=AOArm_grad*(2*pi)/360;

%% CHANNEL MATRIX
% Channel matrix between RIS and BS
Abr=nan(Nr,Lbr); ro_br=zeros(Lbr,Lbr);
Abr_h=nan(Nris,Lbr); % future hermitian
for ll = 1:Lbr
    ro_br(ll,ll)=F*(c/(4*pi*(r_br(ll)+d_br(ll))*(f*10^9)))^(mu/2);
    for b=-B:B
        t=B+b+1;
        Abr(t,ll) = exp(j*(b*omega(AODbr(ll),delta,lambda) + b^(2)*gamma1(AODbr(ll),r_br(ll),delta,lambda)));
    end
    for r=-R:R
        g=R+r+1;
        Abr_h(g,ll) = exp(j*(r*omega(AOAbr(ll),delta, lambda) + r^(2)*gamma1(AOAbr(ll),d_br(ll),delta,lambda)));
    end
end

Hbr = Abr*ro_br*Abr_h';

% Channel matrix between US and RIS
Arm=nan(Nris,Lrm); ro_rm=zeros(Lrm,Lrm);
Arm_f=zeros(Nt,Lrm); % hermitian
for ll = 1:Lrm
    ro_rm(ll,ll)=F*(c/(4*pi*(r_rm(ll)+d_rm(ll))*f))^(mu/2);
    for r=-R:R
        g=R+r+1;
        Arm(g,ll) = exp(j*(r*omega(AODrm(ll),delta,lambda) + r^(2)*gamma1(AODrm(ll),r_rm(ll),delta,lambda)));
    end
    for m=-M:M
        s=M+m+1;
        Arm_f(s,ll) = exp(j*(m*omega(AOArm(ll),delta,lambda) + m^(2)*gamma1(AOArm(ll),d_rm(ll),delta,lambda)));
    end
end

% RIS
RIS=zeros(Nris,Nris);
ampRIS=ones(Nris,Nris);
for r=1:Nris
    RIS(r,r)=ampRIS(r)*exp(j*phaseRIS(r));
end

%% CRLB
% INPUT DATA
eta=[d_rm, r_rm AODrm,AOArm];

% Derivadas de A respecte del angle
[dA_AODrm]=derA_angle('RM_AOD',AODrm_grad, r_rm,Arm,delta,lambda, Nr, Nt,Nris, Lrm);
[dA_AOArm]=derA_angle('RM_AOA',AOArm_grad, d_rm,Arm_f,delta,lambda, Nr, Nt,Nris, Lrm);

% Derivades de A respecte la distancia
[dA_rrm_AOD]=derA_r('RM_AOD', AODrm_grad, r_rm, Arm, delta, lambda, Nr, Nt, Nris, Lbr);
[dA_drm_AOA]=derA_r('RM_AOA', AOArm_grad, d_rm, Arm_f, delta,lambda, Nr, Nt, Nris, Lbr);

% Derivades de RO respecte la distancia
[d_RO]=derRO_r(r_rm, d_rm,mu, c, F, f);

% Expressions de les derivades
dd_rm=Hbr*RIS*(Arm*d_RO*(Arm_f)'+Arm*ro_rm*(dA_drm_AOA)');
dr_rm= Hbr*RIS*(dA_rrm_AOD*ro_rm*(Arm_f)'+Arm*d_RO*(Arm_f)');
dAoDrm= Hbr*RIS*(dA_AODrm*ro_rm*(Arm_f)');
dAoArm= Hbr*RIS*(Arm*ro_rm*(dA_AOArm)');

% passar snr a lineal
SNR_mag=db2pow(SNR);
E4 = 3*(Q/2)^4;
%E4=1;
FIM=nan(length(eta),length(eta));
FIM(1,1)=E4*((SNR_mag)/(Nt))*real(trace(dr_rm'*dr_rm));
FIM(1,2)=E4*((SNR_mag)/(Nt))*real(trace(dr_rm'*dd_rm));
FIM(1,3)=E4*((SNR_mag)/(Nt))*real(trace(dr_rm'*dAoDrm));
FIM(1,4)=E4*((SNR_mag)/(Nt))*real(trace(dr_rm'*dAoArm)); 

FIM(2,1)=E4*((SNR_mag)/(Nt))*real(trace(dd_rm'*dr_rm));
FIM(2,2)=E4*((SNR_mag)/(Nt))*real(trace(dd_rm'*dd_rm));
FIM(2,3)=E4*((SNR_mag)/(Nt))*real(trace(dd_rm'*dAoDrm));
FIM(2,4)=E4*((SNR_mag)/(Nt))*real(trace(dd_rm'*dAoArm)); 

FIM(3,1)=E4*((SNR_mag)/(Nt))*real(trace(dAoDrm'*dr_rm));
FIM(3,2)=E4*((SNR_mag)/(Nt))*real(trace(dAoDrm'*dd_rm));
FIM(3,3)=E4*((SNR_mag)/(Nt))*real(trace(dAoDrm'*dAoDrm));
FIM(3,4)=E4*((SNR_mag)/(Nt))*real(trace(dAoDrm'*dAoArm));

FIM(4,1)=E4*((SNR_mag)/(Nt))*real(trace(dAoArm'*dr_rm)); 
FIM(4,2)=E4*((SNR_mag)/(Nt))*real(trace(dAoArm'*dd_rm)); 
FIM(4,3)=E4*((SNR_mag)/(Nt))*real(trace(dAoArm'*dAoDrm));
FIM(4,4)=E4*((SNR_mag)/(Nt))*real(trace(dAoArm'*dAoArm));

% Transformation matrix
T=nan(length(eta),2);

T(1,1)=-(-xM+xRIS)/(d_rm);
T(1,2)=(-yM+yRIS)/d_rm;

T(2,1)=-(-xM+xRIS)/(d_rm);
T(2,2)=(-yM+yRIS)/d_rm;

T(3,1)=-(yRIS-yM)/(((-xM+xRIS)^2)+(-yM+yRIS)^2);
T(3,2)=(xRIS-xM)/(((-xM+xRIS)^2)+(-yM+yRIS)^2);

T(4,1)=-(yRIS-yM)/(((-xM+xRIS)^2)+(-yM+yRIS)^2);
T(4,2)=(xRIS-xM)/(((-xM+xRIS)^2)+(-yM+yRIS)^2);

FIM_T=T'*FIM*T;
CRLB=((trace(inv(FIM_T))));
end

