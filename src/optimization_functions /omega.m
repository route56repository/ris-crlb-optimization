%% OMEGA
function [om]=omega(angle,delta,lambda)
om=-(2*pi*delta/lambda)*(sin(angle));
end