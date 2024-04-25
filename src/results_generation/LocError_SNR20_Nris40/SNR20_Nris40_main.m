%% Evolution of the RMSE throughout iterations (with NR = 40 and SNR = 20dB) 
% for the three methods and Evolution of the SNR throughout iterations 
% (with NR = 40 and SNR = 20dB) for the three methods.
addpath('/Users/carla.macias/Desktop/ris-crlb-optimization/src/optimization_functions');
warning('off','all')
tic
xBS=0; yBS=0; xRIS=5; yRIS=5; xM_real=10; yM_real=0.5;
Nr=51; Nt= 21; Nris=40;

% SNR definition
Bk=10e6; Tk=290; Boltzman=physconst('Boltzmann');
sigma2=Bk*Tk*Boltzman; M0=64;
SNR=20; SNRlineal=db2pow(SNR);
powerP=SNRlineal*(sigma2*M0);

f = 28; j=sqrt(-1); Q=0.2;
c = physconst('LightSpeed'); lambda = c/(f*10^9);
Lbr=1; Lrm=1; delta=lambda/2; mu=3; B=(Nr-1)/2; M=(Nt-1)/2;
B_vect=-B:B; M_vect=-M:M; R=(Nris-1)/2;

% Iterations parameter
iterations=1;  itMC=2;

%Results vectors definition
MSE_x_vect_CRLB_MC=nan(itMC, iterations); MSE_y_vect_CRLB_MC=nan(itMC, iterations);
MSE_x_vect_SNR_MC=nan(itMC, iterations);  MSE_y_vect_SNR_MC=nan(itMC, iterations);
MSE_x_vect_rand_MC=nan(itMC, iterations); MSE_y_vect_rand_MC=nan(itMC, iterations);

MSE_phases_CRLB=nan(itMC, Nris);
MSE_phases_SNR=nan(itMC, Nris);
MSE_phases_rand=nan(itMC, Nris);

MC_d_rm_estimate_vect_CRLB=[]; MC_r_rm_estimate_vect_CRLB=[];
MC_AODrm_estimate_grad_vect_CRLB=[]; MC_AOArm_estimate_grad_vect_CRLB=[];

RMSEposc_vect_CRLB_MC=nan(itMC, iterations);
RMSEposc_vect_SNR_MC=nan(itMC, iterations);
RMSEposc_vect_rand_MC=nan(itMC, iterations);
MC_iterCRLB=nan(itMC, iterations-1);

MSE_drm_vect_CRLB_MC=nan(itMC, iterations);
MSE_rrm_vect_CRLB_MC=nan(itMC, iterations);
MSE_AODrm_vect_CRLB_MC=nan(itMC, iterations);
MSE_AOArm_vect_CRLB_MC=nan(itMC, iterations);
SNRreception_vect_it_CRLB=nan(itMC, iterations);
SNRreception_vect_it_SNR=nan(itMC, iterations);
SNRreception_vect_it_rand=nan(itMC, iterations);


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

%%
for MC=1:itMC
    % Random paramenters
    F=sqrt(Q/2)*(randn(1)+j*randn(1));
    %phasesRIS_0=2*pi*rand(1,Nris)-pi;
    phasesRIS_0=zeros(1,Nris);

    %% H channel
    % Channel matrix between RIS and BS
    Abr=nan(Nr,Lbr); ro_br=zeros(Lbr,Lbr);
    Abr_h=nan(Nris,Lbr); % future hermitian
    for ll = 1:Lbr
        ro_br(ll,ll)=F*(c/(4*pi*(r_br(ll)+d_br(ll))*(f*10^9)))^(mu/2);
        for b=-B:B
            t=B+b+1;
            Abr(t,ll) = exp(j*(b*omega(AODbr(ll),delta,lambda) + (b^2)*gamma1(AODbr(ll),r_br(ll),delta,lambda)));
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
    Hrm = Arm*ro_rm*Arm_f';


    % X signal
    x = generate_random_prs_matrix(Nt,M0,powerP);

    [SNRreception_vect_SNR,RMSEposc_SNR_vect,MSE_drm_vect_SNR,MSE_rrm_vect_SNR,MSE_AODrm_vect_SNR,MSE_AOArm_vect_SNR, cota_CRLB_vect_SNR, phasesRIS_vect_SNR_02pi,xM_estimat_vect_SNR,yM_estimat_vect_SNR,MSE_x_vect_SNR,MSE_y_vect_SNR, AOArm_estimate_grad_vect_SNR, AODrm_estimate_grad_vect_SNR, r_rm_estimate_vect_SNR,d_rm_estimate_vect_SNR]=SNR_signalY_IT2(Hbr, Hrm, iterations,xRIS,yRIS,xM_real,yM_real,Nr, Nt, Nris, phasesRIS_0,x,M0, sigma2);
    [SNRreception_vect_CRLB,iterCRLB,RMSEposc_CRLB_vect,MSE_drm_vect_CRLB,MSE_rrm_vect_CRLB,MSE_AODrm_vect_CRLB,MSE_AOArm_vect_CRLB, cota_CRLB_vect_CRLB, phasesRIS_vect_CRLB_02pi,xM_estimat_vect_CRLB,yM_estimat_vect_CRLB,MSE_x_vect_CRLB,MSE_y_vect_CRLB, AOArm_estimate_grad_vect_CRLB, AODrm_estimate_grad_vect_CRLB, r_rm_estimate_vect_CRLB,d_rm_estimate_vect_CRLB]=CRLB_signalY_IT_startSNR_noML(Hbr, Hrm, iterations,xRIS,yRIS,xM_real,yM_real,Nr, Nt, Nris, phasesRIS_0, SNR,x,M0,d_rm,r_rm,AODrm_grad,AOArm_grad, sigma2);
    [SNRreception_vect_rand ,RMSEposc_rand_vect,phasesRIS_vect_rand, Ea_x_vect_rand, Ea_y_vect_rand,xM_estimat_vect_rand,yM_estimat_vect_rand,MSE_x_vect_rand,MSE_y_vect_rand, AOArm_estimate_grad_vect_rand, AODrm_estimate_grad_vect_rand, r_rm_estimate_vect_rand]=rand_signalY_IT(Hbr, Hrm, iterations,xRIS,yRIS,xM_real,yM_real,Nr, Nt, Nris, phasesRIS_0, SNR,x,M0, sigma2);

    MC_iterCRLB(MC,:)=iterCRLB;
    MSE_x_vect_CRLB_MC(MC,:)= MSE_x_vect_CRLB;
    MSE_y_vect_CRLB_MC(MC,:)= MSE_y_vect_CRLB;
    RMSEposc_vect_CRLB_MC(MC,:)= RMSEposc_CRLB_vect;
    SNRreception_vect_it_CRLB(MC,:)=SNRreception_vect_CRLB;

    MSE_x_vect_SNR_MC(MC,:)=MSE_x_vect_SNR;
    MSE_y_vect_SNR_MC(MC,:)=MSE_y_vect_SNR;
    RMSEposc_vect_SNR_MC(MC,:)=RMSEposc_SNR_vect;
    SNRreception_vect_it_SNR(MC,:)=SNRreception_vect_SNR;

    MSE_x_vect_rand_MC(MC,:)=MSE_x_vect_rand;
    MSE_y_vect_rand_MC(MC,:)=MSE_y_vect_rand;
    RMSEposc_vect_rand_MC(MC,:)=RMSEposc_rand_vect;
    SNRreception_vect_it_rand(MC,:)=SNRreception_vect_rand;

    MSE_drm_vect_CRLB_MC(MC,:)=MSE_drm_vect_CRLB;
    MSE_rrm_vect_CRLB_MC(MC,:)= MSE_rrm_vect_CRLB;
    MSE_AODrm_vect_CRLB_MC(MC,:)=MSE_AODrm_vect_CRLB;
    MSE_AOArm_vect_CRLB_MC(MC,:)=MSE_AOArm_vect_CRLB;


    %    Needed for bias calculation
    MC_d_rm_estimate_vect_CRLB=[MC_d_rm_estimate_vect_CRLB; d_rm_estimate_vect_CRLB];
    MC_r_rm_estimate_vect_CRLB=[MC_r_rm_estimate_vect_CRLB; r_rm_estimate_vect_CRLB];
    MC_AODrm_estimate_grad_vect_CRLB=[MC_AODrm_estimate_grad_vect_CRLB; AODrm_estimate_grad_vect_CRLB];
    MC_AOArm_estimate_grad_vect_CRLB=[MC_AOArm_estimate_grad_vect_CRLB; AOArm_estimate_grad_vect_CRLB];
end


%% MC
[A,B2]=size(MSE_x_vect_SNR_MC);
MC_MSE_result_x_CRLB=nan(1,B2); MC_MSE_result_y_CRLB=nan(1,B2);
MC_MSE_result_x_SNR=nan(1,B2); MC_MSE_result_y_SNR=nan(1,B2);
MC_MSE_result_x_rand=nan(1,B2); MC_MSE_result_y_rand=nan(1,B2);
MC_MSE_result_drm=nan(1,B2); MC_MSE_result_rrm=nan(1,B2);
MC_MSE_result_AODrm=nan(1,B2); MC_MSE_result_AOArm=nan(1,B2);
MC_RMSE_result_CRLB=nan(1,B2);
MC_RMSE_result_SNR=nan(1,B2);
MC_RMSE_result_rand=nan(1,B2);
MC_SNR_reception_CRLB=nan(1,B2);
MC_SNR_reception_SNR=nan(1,B2);
MC_SNR_reception_rand=nan(1,B2);

MC_result_drm=nan(1,B2); MC_result_rrm=nan(1,B2);
MC_result_AODrm=nan(1,B2); MC_result_AOArm=nan(1,B2);

for comp=1:B2
    MC_MSE_result_x_CRLB(1,comp)=sum(MSE_x_vect_CRLB_MC(:,comp))/A;
    MC_MSE_result_y_CRLB(1,comp)=sum(MSE_y_vect_CRLB_MC(:,comp))/A;
    MC_RMSE_result_CRLB(1,comp)=sum(RMSEposc_vect_CRLB_MC(:,comp))/A;
    MC_SNR_reception_CRLB(1,comp)=sum(SNRreception_vect_it_CRLB(:,comp))/A;

    MC_MSE_result_x_SNR(1,comp)=sum(MSE_x_vect_SNR_MC(:,comp))/A;
    MC_MSE_result_y_SNR(1,comp)=sum(MSE_y_vect_SNR_MC(:,comp))/A;
    MC_RMSE_result_SNR(1,comp)=sum(RMSEposc_vect_SNR_MC(:,comp))/A;
    MC_SNR_reception_SNR(1,comp)=sum(SNRreception_vect_it_SNR(:,comp))/A;

    MC_MSE_result_x_rand(1,comp)=sum(MSE_x_vect_rand_MC(:,comp))/A;
    MC_MSE_result_y_rand(1,comp)=sum(MSE_y_vect_rand_MC(:,comp))/A;
    MC_RMSE_result_rand(1,comp)=sum(RMSEposc_vect_rand_MC(:,comp))/A;
    MC_SNR_reception_rand(1,comp)=sum(SNRreception_vect_it_rand(:,comp))/A;

    MC_MSE_result_drm(1,comp)=sum(MSE_drm_vect_CRLB_MC(:,comp))/A;
    MC_MSE_result_rrm(1,comp)=sum(MSE_rrm_vect_CRLB_MC(:,comp))/A;
    MC_MSE_result_AODrm(1,comp)=sum(MSE_AODrm_vect_CRLB_MC(:,comp))/A;
    MC_MSE_result_AOArm(1,comp)=sum(MSE_AOArm_vect_CRLB_MC(:,comp))/A;

    MC_result_drm(1,comp)=sum(MC_d_rm_estimate_vect_CRLB(:,comp))/A;
    MC_result_rrm(1,comp)=sum(MC_r_rm_estimate_vect_CRLB(:,comp))/A;
    MC_result_AODrm(1,comp)=sum(MC_AODrm_estimate_grad_vect_CRLB(:,comp))/A;
    MC_result_AOArm(1,comp)=sum(MC_AOArm_estimate_grad_vect_CRLB(:,comp))/A;

end

%%
[A3,B3]=size(MC_iterCRLB);
MC_iterCRLB_result=nan(1,B3);
for comp3=1:B3
    MC_iterCRLB_result(1,comp3)=sum(MC_iterCRLB(:,comp3))/A3;
end

RMSE2_CRLB=sqrt(MC_MSE_result_x_CRLB+MC_MSE_result_y_CRLB);
RMSE2_SNR=sqrt(MC_MSE_result_x_SNR+MC_MSE_result_y_SNR);
RMSE2_rand=sqrt(MC_MSE_result_x_rand+MC_MSE_result_y_rand);
%%
% CRLB 
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

func=@(phaseRIS)CRLBfunction(xBS,yBS,xRIS,yRIS,xM_real,yM_real,phaseRIS,Nris,SNR, Q);
options = optimset('TolX',1e-6,'TolFun',1e-6, 'MaxFunEvals', 1000000*Nris, 'MaxIter', 1000000*Nris);
[fases,fval]=fminsearch(func,phasesRIS_SNR);
MC_cotaCRLB=ones(1,iterations)*sqrt(fval);


%%

%Bias calculation, miro que % del MSE es bias, como mi estimador es no sesgado tendría que ser casi 0
bias_vect_d_rm=(d_rm-MC_result_drm).^2;
bias_vect_r_rm=(r_rm-MC_result_rrm).^2;
bias_vect_AODrm=(AODrm_grad-MC_result_AODrm).^2;
bias_vect_AOArm=(AOArm_grad-MC_result_AOArm).^2;

comp_biasMSEdrm=100*bias_vect_d_rm(end)/MC_MSE_result_drm(end);
comp_biasMSErrm=100*bias_vect_r_rm(end)/MC_MSE_result_rrm(end);
comp_biasMSEAODrm=100*bias_vect_AODrm(end)/MC_MSE_result_AODrm(end);
comp_biasMSEAOArm=100*bias_vect_AOArm(end)/MC_MSE_result_AOArm(end);


%%
fig1=figure(1);
semilogy(MC_RMSE_result_CRLB,'-ob','LineWidth',1)
hold on
semilogy(MC_RMSE_result_SNR,'-og','LineWidth',1)
hold on
semilogy(MC_RMSE_result_rand,'-or','LineWidth',1)
hold on
semilogy(1:iterations,MC_cotaCRLB,'--r','LineWidth',1)
grid on
legend('Min. CRLB phases method','Max. SNR phases method','Random Phases','CRLB','CRLB evolution','Location','southwest')
xlabel('Iterations')
ylabel('RMSE')
hold off

fig2=figure(2);
plot(1:iterations, MC_SNR_reception_CRLB,'LineWidth',1)
hold on
plot(1:iterations, MC_SNR_reception_SNR,'LineWidth',1)
hold on
plot(1:iterations, MC_SNR_reception_rand,'LineWidth',1)
grid on
xlabel('Iterations')
ylabel('SNR')
legend('CRLB','SNR','RAND')
hold off
toc
%save('SNR20_Nr40')

