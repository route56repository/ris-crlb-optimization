
%% RMSE as a function of the RIS size for the two studied methods with SNR = 20.

clear all
% Vector definitions
SNR_vect=-20:5:20;
error_CRLB=nan(1,length(SNR_vect));
error_SNR=nan(1,length(SNR_vect));
load('LocErrDiffSNR_Nris40.mat')
close all

% RMSE for SNR=-20
MC_RMSE_result_CRLB=mean(RMSEposcCRLB_m20);
error_CRLB(1)=MC_RMSE_result_CRLB(end);
MC_RMSE_result_SNR=mean(RMSEposcSNRm20);
error_SNR(1)=MC_RMSE_result_SNR(end);

% RMSE for SNR=-15
MC_RMSE_result_CRLB=mean(RMSEposcCRLB_m15);
error_CRLB(2)=MC_RMSE_result_CRLB(end);
MC_RMSE_result_SNR=mean(RMSEposcSNR_m15);
error_SNR(2)=MC_RMSE_result_SNR(end);

% RMSE for SNR=-10
MC_RMSE_result_CRLB=mean(RMSEposcCRLB_m10);
error_CRLB(3)=MC_RMSE_result_CRLB(end);
MC_RMSE_result_SNR=mean(RMSEposcSNR_m10);
error_SNR(3)=MC_RMSE_result_SNR(end);

% RMSE for SNR=-5
MC_RMSE_result_CRLB=mean(RMSEposcCRLB_m5);
error_CRLB(4)=MC_RMSE_result_CRLB(end);
MC_RMSE_result_SNR=mean(RMSEposcSNR_m5);
error_SNR(4)=MC_RMSE_result_SNR(end);

% RMSE for SNR=0
MC_RMSE_result_CRLB=mean(RMSEposcCRLB_0);
error_CRLB(5)=MC_RMSE_result_CRLB(end);
MC_RMSE_result_SNR=mean(RMSEposcSNR_0);
error_SNR(5)=MC_RMSE_result_SNR(end);

% RMSE for SNR=5
MC_RMSE_result_CRLB=mean(RMSEposcCRLB_5);
error_CRLB(6)=MC_RMSE_result_CRLB(end);
MC_RMSE_result_SNR=mean(RMSEposcSNR_5);
error_SNR(6)=MC_RMSE_result_SNR(end);

% RMSE for SNR=10
MC_RMSE_result_CRLB=mean(RMSEposcCRLB_10);
error_CRLB(7)=MC_RMSE_result_CRLB(end);
MC_RMSE_result_SNR=mean(RMSEposcSNR_10);
error_SNR(7)=MC_RMSE_result_SNR(end);

% RMSE for SNR=15
MC_RMSE_result_CRLB=mean(RMSEposcCRLB_15);
error_CRLB(8)=MC_RMSE_result_CRLB(end);
MC_RMSE_result_SNR=mean(RMSEposcSNR_15);
error_SNR(8)=MC_RMSE_result_SNR(end);

% RMSE for SNR=20
MC_RMSE_result_CRLB=mean(RMSEposcCRLB_20);
error_CRLB(9)=MC_RMSE_result_CRLB(end);
MC_RMSE_result_SNR=mean(RMSEposcSNR_20);
error_SNR(9)=MC_RMSE_result_SNR(end);

% Plot result
fig1=figure(1)
SNR_vect=-20:5:20;
semilogy(SNR_vect, error_CRLB,'-ob','LineWidth',1.2)
hold on
semilogy(SNR_vect, error_SNR,'-o','LineWidth',1.2)
grid on
xlabel('SNR')
ylabel('RMSE [m]')
legend('Min. CRLB phases method','Max. SNR phases method')
hold off




