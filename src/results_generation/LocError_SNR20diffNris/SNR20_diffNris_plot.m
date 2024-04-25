
%% RMSE as a function of the RIS size for the two studied methods with SNR = 20.

clear all
% Vector definitions
Nris_vect=[10 20 30 40 50 60 70 80 90 100];
error_CRLB=nan(1,length(Nris_vect));
error_SNR=nan(1,length(Nris_vect));
load('LocErrDiffNris_10to100.mat')
close all

% RMSE for Nris=10
MC_RMSE_result_CRLB=mean(RMSEposcCRLB_10);
error_CRLB(1)=MC_RMSE_result_CRLB(end);
MC_RMSE_result_SNR=mean(RMSEposcSNR_10);
error_SNR(1)=MC_RMSE_result_SNR(end);

% RMSE for Nris=20
MC_RMSE_result_CRLB=mean(RMSEposcCRLB_20);
error_CRLB(2)=MC_RMSE_result_CRLB(end);
MC_RMSE_result_SNR=mean(RMSEposcSNR_20);
error_SNR(2)=MC_RMSE_result_SNR(end);

% RMSE for Nris=30
MC_RMSE_result_CRLB=mean(RMSEposcCRLB_30);
error_CRLB(3)=MC_RMSE_result_CRLB(end);
MC_RMSE_result_SNR=mean(RMSEposcSNR_30);
error_SNR(3)=MC_RMSE_result_SNR(end);

% RMSE for Nris=40
MC_RMSE_result_CRLB=mean(RMSEposcCRLB_40);
error_CRLB(4)=MC_RMSE_result_CRLB(end);
MC_RMSE_result_SNR=mean(RMSEposcSNR_40);
error_SNR(4)=MC_RMSE_result_SNR(end);

% RMSE for Nris=50
MC_RMSE_result_CRLB=mean(RMSEposcCRLB_50);
error_CRLB(5)=MC_RMSE_result_CRLB(end);
MC_RMSE_result_SNR=mean(RMSEposcSNR_50);
error_SNR(5)=MC_RMSE_result_SNR(end);

% RMSE for Nris=60
MC_RMSE_result_CRLB=mean(RMSEposcCRLB_60);
error_CRLB(6)=MC_RMSE_result_CRLB(end);
MC_RMSE_result_SNR=mean(RMSEposcSNR_60);
error_SNR(6)=MC_RMSE_result_SNR(end);

% RMSE for Nris=70
MC_RMSE_result_CRLB=mean(RMSEposcCRLB_70);
error_CRLB(7)=MC_RMSE_result_CRLB(end);
MC_RMSE_result_SNR=mean(RMSEposcSNR_70);
error_SNR(7)=MC_RMSE_result_SNR(end);

% RMSE for Nris=80
MC_RMSE_result_CRLB=mean(RMSEposcCRLB_80);
error_CRLB(8)=MC_RMSE_result_CRLB(end);
MC_RMSE_result_SNR=mean(RMSEposcSNR_80);
error_SNR(8)=MC_RMSE_result_SNR(end);

% RMSE for Nris=90
MC_RMSE_result_CRLB=mean(RMSEposcCRLB_90);
error_CRLB(9)=MC_RMSE_result_CRLB(end);
MC_RMSE_result_SNR=mean(RMSEposcSNR_90);
error_SNR(9)=MC_RMSE_result_SNR(end);

% RMSE for Nris=100
MC_RMSE_result_CRLB=mean(RMSEposcCRLB_100);
error_CRLB(100)=MC_RMSE_result_CRLB(end);
MC_RMSE_result_SNR=mean(RMSEposcSNR_100);
error_SNR(100)=MC_RMSE_result_SNR(end);

% Plot result
fig1=figure(1)
Nris_vect=[10 20 30 40 50 60 70 80 90 100];
semilogy(Nris_vect, error_CRLB,'-ob','LineWidth',1.2)
hold on
semilogy(Nris_vect, error_SNR,'-o','LineWidth',1.2)
grid on
xlabel('N^{R}')
ylabel('RMSE [m]')
legend('Min. CRLB phases method','Max. SNR phases method')
hold off




