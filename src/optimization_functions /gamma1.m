function [gam]=gamma1(angle,dist,delta,lambda)
gam=(pi*delta^(2)/(lambda*dist))*(cos(angle))^2;
end
