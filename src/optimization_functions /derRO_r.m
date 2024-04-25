
function [d_RO]=derRO_r(r_rm, d_rm,mu, c, F, f)
 d_RO=-F*(mu*c^(mu/2))/(2^(mu/2)*pi^(mu/2)*f^(mu/2)*(d_rm+r_rm)^((mu+2)/2));
end