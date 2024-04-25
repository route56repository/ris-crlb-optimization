function [signal_noise,SNRreception_dB]=awgnoise(signal,sigma2)

j=sqrt(-1);
[numSamples,col1]=size(signal);

% % % % % Calculate signal power
% % % % signalPower = mean(abs(signal).^2);
% % % % 
% % % % %noisePower = signalPower / desiredSNR;
% % % % noisePower=power/

% Generate Gaussian noise with the calculated noise power
noise = sqrt(sigma2/2) * (randn(1, col1) + j * randn(1, col1));

 % Add noise to the complex signal
 signal_noise= signal + noise;
 %comp=pow2db((signal*signal')/abs((signal_noise*signal_noise'-signal*signal')));
% %comp=pow2db((signal'*signal)/abs((signal_noise'*signal_noise-signal'*signal)));

signalPower = mean(abs(signal).^2);
noisePower = mean(abs(noise).^2);
SNRreception_dB=pow2db(signalPower/noisePower);
end