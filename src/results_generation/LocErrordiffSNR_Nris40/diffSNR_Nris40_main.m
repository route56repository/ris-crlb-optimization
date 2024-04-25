%% SIMULATION USING THE Y SIGNAL
clear all
addpath('/Users/carla.macias/Desktop/ris-crlb-optimization/src/optimization_functions');
warning('off','all')
xBS=0; yBS=0; xRIS=5; yRIS=5; xM_real=10; yM_real=0.5;
Nr=51; Nt= 21;

% SNR definition
Bk=10e6; Tk=290; Boltzman=physconst('Boltzmann');
sigma2=Bk*Tk*Boltzman; M0=64;
 SNRlineal=db2pow(SNR);
powerP=SNRlineal*(sigma2*M0);

Nris=40;
SNR_vect=-20:5:20;
SNR_RMSE_CRLB=nan(1,length(SNR_vect));
SNR_RMSE_SNR=nan(1,length(SNR_vect));
w=1;
for SNR=SNR_vect
    f = 28; j=sqrt(-1); Q=0.2;
    c = physconst('LightSpeed'); lambda = c/(f*10^9);
    Lbr=1; Lrm=1; delta=lambda/2; mu=3; B=(Nr-1)/2; M=(Nt-1)/2;
    B_vect=-B:B; M_vect=-M:M; R=(Nris-1)/2;

    % Iterations parameter
    iterations=10;  itMC=1000;

    %Results vectors definition
    MSE_x_vect_SNR_MC=nan(itMC, iterations); 
    RMSEposc_vect_CRLB_MC=nan(itMC, iterations);
    RMSEposc_vect_SNR_MC=nan(itMC, iterations);

 
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
   
    for MC=1:itMC
        F=sqrt(Q/2)*(randn(1)+j*randn(1));
        phasesRIS_0=zeros(1,Nris);

        % H channel
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

        x = generate_random_prs_matrix(Nt,M0,powerP);
        [SNRreception_vect_SNR,RMSEposc_SNR_vect,MSE_drm_vect_SNR,MSE_rrm_vect_SNR,MSE_AODrm_vect_SNR,MSE_AOArm_vect_SNR, cota_CRLB_vect_SNR, phasesRIS_vect_SNR_02pi,xM_estimat_vect_SNR,yM_estimat_vect_SNR,MSE_x_vect_SNR,MSE_y_vect_SNR, AOArm_estimate_grad_vect_SNR, AODrm_estimate_grad_vect_SNR, r_rm_estimate_vect_SNR,d_rm_estimate_vect_SNR]=SNR_signalY_IT2(Hbr, Hrm, iterations,xRIS,yRIS,xM_real,yM_real,Nr, Nt, Nris, phasesRIS_0,x,M0, sigma2);
        [SNRreception_vect_CRLB,iterCRLB,RMSEposc_CRLB_vect,MSE_drm_vect_CRLB,MSE_rrm_vect_CRLB,MSE_AODrm_vect_CRLB,MSE_AOArm_vect_CRLB, cota_CRLB_vect_CRLB, phasesRIS_vect_CRLB_02pi,xM_estimat_vect_CRLB,yM_estimat_vect_CRLB,MSE_x_vect_CRLB,MSE_y_vect_CRLB, AOArm_estimate_grad_vect_CRLB, AODrm_estimate_grad_vect_CRLB, r_rm_estimate_vect_CRLB,d_rm_estimate_vect_CRLB]=CRLB_signalY_IT_startSNR_noML(Hbr, Hrm, iterations,xRIS,yRIS,xM_real,yM_real,Nr, Nt, Nris, phasesRIS_0, SNR,x,M0,d_rm,r_rm,AODrm_grad,AOArm_grad, sigma2);
        RMSEposc_vect_CRLB_MC(MC,:)= RMSEposc_CRLB_vect;
        RMSEposc_vect_SNR_MC(MC,:)=RMSEposc_SNR_vect;
    end


    % MC
    MC_RMSE_result_CRLB=mean(RMSEposc_vect_CRLB_MC);
    MC_RMSE_result_SNR=mean(RMSEposc_vect_SNR_MC);


    if SNR==-20
        RMSEposcCRLB_m20=RMSEposc_vect_CRLB_MC;
        RMSEposcSNR_m20=RMSEposc_vect_SNR_MC;
    elseif SNR==-15
        RMSEposcCRLB_m15=RMSEposc_vect_CRLB_MC;
        RMSEposcSNR_m15=RMSEposc_vect_SNR_MC;
    elseif SNR==-10
        RMSEposcCRLB_m10=RMSEposc_vect_CRLB_MC;
        RMSEposcSNR_m10=RMSEposc_vect_SNR_MC;
    elseif SNR==-5
        RMSEposcCRLB_m5=RMSEposc_vect_CRLB_MC;
        RMSEposcSNR_m5=RMSEposc_vect_SNR_MC;
    elseif SNR==0
        RMSEposcCRLB_0=RMSEposc_vect_CRLB_MC;
        RMSEposcSNR_0=RMSEposc_vect_SNR_MC;
    elseif SNR==5
        RMSEposcCRLB_5=RMSEposc_vect_CRLB_MC;
        RMSEposcSNR_5=RMSEposc_vect_SNR_MC;
    elseif SNR==10
        RMSEposcCRLB_10=RMSEposc_vect_CRLB_MC;
        RMSEposcSNR_10=RMSEposc_vect_SNR_MC;
    elseif SNR==15
        RMSEposcCRLB_15=RMSEposc_vect_CRLB_MC;
        RMSEposcSNR_15=RMSEposc_vect_SNR_MC;
    elseif SNR==20
        RMSEposcCRLB_20=RMSEposc_vect_CRLB_MC;
        RMSEposcSNR_20=RMSEposc_vect_SNR_MC;
    end

    SNR_RMSE_CRLB(w)=MC_RMSE_result_CRLB(end);
    SNR_RMSE_SNR(w)=MC_RMSE_result_SNR(end);
    w=w+1;
end

fig1=figure(1);
semilogy(SNR_vect, SNR_RMSE_CRLB,'-ob','LineWidth',1)
hold on
semilogy(SNR_vect, SNR_RMSE_SNR,'-og','LineWidth',1)
hold on
grid on
legend('Min. CRLB phases method','Max. SNR phases method','Location','southwest')
xlabel('SNR [dB]')
ylabel('RMSE [m]')
hold off

save('LocErrDiffSNR_Nris40')

