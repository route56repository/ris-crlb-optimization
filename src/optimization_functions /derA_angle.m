% % % % % % % % function [dA_angle,T]=derA_angle(type,angle, dist,matriu,delta,lambda, Nr, Nt,Nris)
% % % % % % % % j = sqrt(-1);
% % % % % % % % B=(Nr-1)/2; M=(Nt-1)/2; R=(Nris-1)/2;
% % % % % % % % if type=='BR_AOD'
% % % % % % % %     T=nan(Nr,1);
% % % % % % % %     for b=-B:B
% % % % % % % %         t=B+b+1;
% % % % % % % %         T(t)=-((j*b*2*pi*delta)/(lambda))*cosd(angle)+j*b^(2)*((pi*delta^(2))/(lambda*dist))*2*cosd(angle)*(-sind(angle));
% % % % % % % %     end
% % % % % % % % elseif type=='BR_AOA'
% % % % % % % %     T=nan(Nris,1);
% % % % % % % %     for r=-R:R
% % % % % % % %         t=R+r+1;
% % % % % % % %         T(t)=-((j*r*2*pi*delta)/(lambda))*cosd(angle)+j*r^(2)*((pi*delta^(2))/(lambda*dist))*2*cosd(angle)*(-sind(angle));
% % % % % % % %     end
% % % % % % % % elseif type=='RM_AOD'
% % % % % % % %     T=nan(Nris,1);
% % % % % % % %     for r=-R:R
% % % % % % % %         t=R+r+1;
% % % % % % % %         T(t)=-((j*r*2*pi*delta)/(lambda))*cosd(angle)+j*r^(2)*((pi*delta^(2))/(lambda*dist))*2*cosd(angle)*(-sind(angle));
% % % % % % % %     end
% % % % % % % % 
% % % % % % % % elseif type=='RM_AOA'
% % % % % % % %     T=nan(Nt,1);
% % % % % % % %     for m=-M:M
% % % % % % % %         t=M+m+1;
% % % % % % % %         T(t)=-((j*m*2*pi*delta)/(lambda))*cosd(angle)+j*m^(2)*((pi*delta^(2))/(lambda*dist))*2*cosd(angle)*(-sind(angle));
% % % % % % % %     end
% % % % % % % % end
% % % % % % % % dA_angle=T.*matriu;
% % % % % % % % end

% MODIFICAT
function [dA_angle,T1]=derA_angle(type,angle, dist,matriu,delta,lambda, Nr, Nt,Nris,Lrm)
j = sqrt(-1);
M=(Nt-1)/2; R=(Nris-1)/2;

if type=='RM_AOD'
    T1=nan(Nris,Lrm);
    for r=-R:R
        t=R+r+1;
        T1(t)=j*(-(2*pi*r*delta*cosd(angle))/(lambda) - (delta^(2)*pi*r^(2)*2*cosd(angle)*sind(angle))/(lambda*dist));
    end

elseif type=='RM_AOA'
    T1=nan(Nt,Lrm);
    for m=-M:M
        t=M+m+1;
        T1(t)=j*(-(2*pi*m*delta*cosd(angle))/(lambda) - (delta^(2)*pi*m^(2)*2*cosd(angle)*sind(angle))/(lambda*dist));
    end
end
dA_angle=T1.*matriu;
end