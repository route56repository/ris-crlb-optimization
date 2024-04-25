
function [dA_r, G]=derA_r(type,angle, dist,matriu,delta,lambda, Nr, Nt, Nris, Lbr)
j = sqrt(-1);
B=(Nr-1)/2; M=(Nt-1)/2; R=(Nris-1)/2;
if type=='RM_AOD'
    G=nan(Nris,Lbr);
    for r=-R:R
        t=R+r+1;
        G(t)=-((j*r^(2)*pi*(cosd(angle))^(2)*delta^(2))/(lambda*dist^2));
    end

elseif type=='RM_AOA'
    G=nan(Nt,Lbr);
    for m=-M:M
        t=M+m+1;
        G(t)=-((j*m^(2)*pi*(cosd(angle))^(2)*delta^(2))/(lambda*dist^2));
    end
end
dA_r=G.*matriu;
end