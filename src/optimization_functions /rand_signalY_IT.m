%% RIS WITH RANDOM PHASE
function  [SNRreception_vect_rand,RMSEposc_rand_vect, phasesRIS_vect_rand, Ea_x_vect_rand, Ea_y_vect_rand,xM_estimat_vect_rand,yM_estimat_vect_rand,MSE_x_vect_rand,MSE_y_vect_rand, AOArm_estimate_grad_vect_rand, AODrm_estimate_grad_vect_rand, r_rm_estimate_vect_rand]=rand_signalY_IT(Hbr, Hrm, iterations,xRIS,yRIS,xM_real,yM_real,Nr, Nt, Nris, phasesRIS_0, SNR,x,M0, sigma2)
B=(Nr-1)/2; M=(Nt-1)/2;
f = 28; c = physconst('LightSpeed'); lambda = c/(f*10^9);
Lbr=1; Lrm=1; delta=lambda/2; j=sqrt(-1);
N1=1000;


MSE_x_vect_rand=[]; MSE_y_vect_rand=[];
xM_estimat_vect_rand=[]; yM_estimat_vect_rand=[];
Ea_x_vect_rand=[]; Ea_y_vect_rand=[];

AOArm_estimate_grad_vect_rand=[]; AODrm_estimate_grad_vect_rand=[]; r_rm_estimate_vect_rand=[];
phasesRIS_vect_rand=phasesRIS_0;
RMSEposc_rand_vect=[];


comp_vect=[];
SNRreception_vect_rand=nan(1,iterations);

for it=1:iterations
    phasesRIS=zeros(1,Nris);
    
     %RIS
    RIS=zeros(Nris,Nris);
    ampRIS=ones(Nris,Nris);
    for r=1:Nris
        RIS(r,r)=ampRIS(r)*exp(j*phasesRIS(r));
    end

    Hf_Y=Hbr*RIS*Hrm;
    Y2=Hf_Y*x;

    Y=nan(Nr,M0);
    for u=1:Nr
        [Y(u,:),SNRreception_dB]=awgnoise(Y2(u,:),sigma2);
        comp_vect=[comp_vect (SNRreception_dB)];
    end
    SNRreception_vect_rand(1,it)=[mean(comp_vect)];

    Happrox=Y*pinv(x);

    %% ESTIMATION OF AODbr
    % V MATRIX
    V=nan(B,M);
    for b=1:B
        t=B+b+1; pt=Nr-t+1;
        for m=1:M
            s=M+m+1; ps=Nt-s+1;
            V(b,m)=Happrox(t,s)*conj(Happrox(pt,ps));
        end
    end

    % S MATRIX
    S1=nan(B,N1);
    for b=1:B
        for n=1:N1
            ang=(2*pi*n-pi*(N1+1))/(N1-1);
            S1(b,n)=exp(2*j*b*omega(ang,delta,lambda));
        end
    end
    Phi1=S1; K=Lbr;
    [q,s] = size(V);
    C1=nan(N1,M); alpha=0.5;
    for k = 1:s
        v=V(:,k);
        C = OMPa(Phi1,v,K,alpha);
        C1(:,k) = C;
    end
    [row,col]=find(C1);
    AODbr_estimate=((2*pi*row-pi*(N1+1))/(N1-1));

    AODbr_estimate_grad1=AODbr_estimate*360/(2*pi);
    AODbr_estimate_grad=AODbr_estimate_grad1(1:Lbr);

    AODbr_estimate_rad=AODbr_estimate_grad*(2*pi)/360;

    %% ESTIMATION OF AOAbr (ho puc fer aixi prk no considero scatters, es a dir, nomes un path)
    AOAbr_estimate_grad=nan(1,Lbr);
    for i=1:length(AODbr_estimate_grad)
        AOAbr_estimate_grad(i)=180+AODbr_estimate_grad(i);
    end
    AOAbr_estimate_rad=AOAbr_estimate_grad*(2*pi)/360;

    %% ESTIMATION OF r_br
    % S3 MATRIX
    Abr3=nan(Nr,N1); r_br_estimate=nan(1,Lbr);
    D=delta*(Nris); d=(2*D^2)/lambda; dist=linspace(0,d,N1);
    for i=1:Lbr
        for ll = 1:N1
            for b=-B:B
                t=B+b+1;
                Abr3(t,ll) = exp(j*(b*omega(AODbr_estimate_rad(i),delta,lambda)+b^(2)*gamma1(AODbr_estimate_rad(i),dist(ll),delta,lambda)));
            end
        end

        Phi3=Abr3; K=Lbr; % sparsity
        [q3,s3] = size(Y); alpha=0.5;
        C3=nan(N1,Nt);
        for k = 1:s3
            v=Y(:,k);
            C = OMPa(Phi3,v,K,alpha);
            C3(:,k) = C;
        end
        [row,col]=find(C3);
        r_br_estimate(i)=dist(row(1));
    end
    d_br_estimate=r_br_estimate;
    %% ESTIMATION OF AOArm
    V2=V';
    % S2 MATRIX
    S2=nan(M,N1);
    for m=1:M
        for n=1:N1
            ang=(2*pi*n-pi*(N1+1))/(N1-1);
            S2(m,n)=exp(-2*j*m*omega(ang,delta,lambda));
        end
    end
    Phi2=S2; K=Lrm; % sparsity
    [q2,s2] = size(V2);
    C2=nan(N1,M); alpha=0.5;
    for k = 1:s2
        v=V2(:,k);
        [C] = orthogonalMatchingPursuit(Phi2, v, K);
        C2(:,k) = C;
    end
    [row,col]=find(C2);
    AOArm_estimate=((2*pi*row(1)-pi*(N1+1))/(N1-1));
    AOArm_estimate_grad1=AOArm_estimate*360/(2*pi);

    [AOArm_estimate_grad]=thirdquadrant(AOArm_estimate_grad1);
    AOArm_estimate_rad=AOArm_estimate_grad*(2*pi)/360;
   
    %% ESTIMATION OF AODrm
    AODrm_estimate_grad=AOArm_estimate_grad-180;
    AODrm_estimate_rad=AODrm_estimate_grad*(2*pi)/360;

    %% ESTIMATION OF d_rm
    V4=Y';
    %S4 MATRIX
    D=delta*(Nris);
    N1=1000;
    Arm4=nan(Nt,N1); d=(2*D^2)/lambda; dist=linspace(0,d,N1);
   
    for  ll= 1:N1
        for m=-M:M
            s4=M+m+1;
            Arm4(s4,ll) = exp(j*(m*omega(AOArm_estimate_rad,delta,lambda) + m^(2)*gamma1(AOArm_estimate_rad,dist(ll),delta,lambda)));
        end
    end
    Phi4=x'*Arm4; 
    K=1; % sparsity
    [q4,s4] = size(V4); alpha=0.6;
    C4=nan(N1,Nr);
    for k = 1:s4
        v=V4(:,k);
        C = OMPa(Phi4,v,K,alpha);
        C4(:,k) = C;
    end

    [row,col]=find(C4);
    d_rm_estimate=dist((row(1)));
    r_rm_estimate=d_rm_estimate;   

    %% POSITION ESTIMATION   
    xM_estimat=xRIS+d_rm_estimate*abs(sind(AODrm_estimate_grad)); %abs prq com esta al primer quadrant el sin sempre es >0
    yM_estimat=yRIS-d_rm_estimate*abs(cosd(AODrm_estimate_grad));
    
    % RIS phasesauan
    MSE_x=abs(xM_real-xM_estimat)^2;
    MSE_y=abs(yM_real-yM_estimat)^2;
    MSE_x_vect_rand=[MSE_x_vect_rand MSE_x];
    MSE_y_vect_rand=[MSE_y_vect_rand MSE_y];

    RMSEposc=sqrt(MSE_x+MSE_y);
    RMSEposc_rand_vect=[RMSEposc_rand_vect RMSEposc];

    Ea_x=abs(xM_real-xM_estimat);
    Ea_y=abs(yM_real-yM_estimat);

    Ea_x_vect_rand=[Ea_x_vect_rand Ea_x];
    Ea_y_vect_rand=[Ea_y_vect_rand Ea_y];

    xM_estimat_vect_rand=[xM_estimat_vect_rand xM_estimat];
    yM_estimat_vect_rand=[yM_estimat_vect_rand yM_estimat];
end
end