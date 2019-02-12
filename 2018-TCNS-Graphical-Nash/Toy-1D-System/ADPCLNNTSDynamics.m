function Zdot=ADPCLNNTSDynamics(t,Z,Zdim,A,A0,Nset,P,L,M,e_histC,...
    x_histC,x0_histC,xDijC,xDi0C,n,m,RC,QC,etac1,etac2,etaa1,...
    etaa2,beta,nu,GammaBar,kmat,GammaThetaC,kTheta,fs,F,g,M_theta,N)
fprintf('%f',t);
global EC;
global UC;
global MUC;
x0te = exp(-0.1*t);
[SC,ASC,A0SC,~,NmC] = ADPCLNNTSGenerateSubgraphs(A,A0,Nset);
ZC=mat2cell(Z,Zdim);
xC=mat2cell(ZC{1},n*ones(N,1));
xHC=mat2cell(ZC{2},n*ones(N,1));
vecthetaHC=mat2cell(ZC{3},n*P,1);
thetaHC=cell(N,1);
WcHC=mat2cell(ZC{4},L);
WaHC=mat2cell(ZC{5},L);
GammaC=mat2cell(ZC{6},L.^2);
% E_i_kC: Cell array containing extrapolation points E_i_k
% x_i_kC: Cell array containing extrapolation points x_i_k
eC=cell(N,1);
for i=1:N
    Si = SC{i};
    ASi = ASC{i};
%     ASiC = ASC(Si);
    A0Si = A0SC{i};
%     if A0Si(1,1)==0
%         x0=zeros(n,1);
%     else
        x0=x0te;
%     end
%     A0SiC = A0SC(Si);
    xSiC = xC(Si);
    Nmi = NmC{i};
%     NmSiC = NmC(Si);
%     mSi = m(Si);
    xDijSiC = xDijC(Si,Si);
    xDi0SiC = xDi0C(Si);
    xDi0 = xDi0SiC{1};
%     thetaH_i = thetaHC{i};
%     thetaHSiC = thetaHC(Si);
    eC{i}=ADPCLNNTSe(Si,xSiC,ASi,A0Si,xDijSiC,xDi0,Nmi,x0);
    thetaHC{i}=reshape(vecthetaHC{i},P(i),n);
end
% EC=cell(length(N),1);
WcH_DotC=cell(N,1);
WaH_DotC=cell(N,1);
xH_DotC=cell(N,1);
vecthetaH_DotC=cell(N,1);
SIDConLearC=cell(N,1);
SIDRankC=cell(N,1);
x_DotC=cell(N,1);
Gamma_DotC=cell(N,1);
uC=cell(N,1);
muC=cell(N,1);
for i=1:N
    Si = SC{i};
    x_i = xC{i};
    e_i = eC{i};
    ASi = ASC{i};
%     ASiC = ASC(Si);
    A0Si = A0SC{i};
%     if A0Si(1,1)==0
%         x0=zeros(n,1);
%     else
        x0=x0te;
%     end
%     A0SiC = A0SC(Si);
    xSiC = xC(Si);
    eSiC = eC(Si);
%     xDi0 = xDi0SiC{1};
    RSiC = RC(Si);
    Qi = QC{i};
    R_i = RC{i};

%     Nmi = NmC{i};
    NmSiC = NmC(Si);
    mSi = m(Si);
    xDijSiC = xDijC(Si,Si);
    xDi0SiC = xDi0C(Si);
%     xDi0 = xD0SiC{i};
    thetaH_i = thetaHC{i};
    xH_i = xHC{i};
    thetaHSiC = thetaHC(Si);
    WcH_i=WcHC{i};
    WaH_i=WaHC{i};
    WaHSiC=WaHC(Si);
    Gamma_i=reshape(GammaC{i},L(i),L(i));
%     ESi=[cell2mat(e(Si));xSiC{1}];
%     EC{i}=ESi;
    
    eta_c1i=etac1(i);
    eta_c2i=etac2(i);
    eta_a1i=etaa1(i);
    eta_a2i=etaa2(i);
    beta_i=beta(i);
    nu_i=nu(i);
    GammaBar_i=GammaBar(i);

    FHi=ADPCLNNTSFHi(xSiC,thetaHSiC,x0,Si,NmSiC,ASi,A0Si,xDijSiC,xDi0SiC,t);
    muHSiC = ADPCLNNTSmuSi(Si,xSiC,eSiC,mSi,NmSiC,ASi,A0Si,xDijSiC,RSiC,WaHSiC,Nset,n);
    omega_i = ADPCLNNTSOmegai(Si,xSiC,eSiC,thetaHSiC,WaHSiC,x0,mSi,NmSiC,ASi,A0Si,xDijSiC,xDi0SiC,RSiC,Nset,n,t);
    GSigmai = ADPCLNNTSGSigmai(Si,xSiC,eSiC,mSi,NmSiC,ASi,A0Si,xDijSiC,Nset,n);
    rho_i = 1+nu_i*omega_i'*Gamma_i*omega_i;
    L_gi = ADPCLNNTSLgi(Si,xSiC,xDijSiC,ASi,A0Si,mSi);
%     muHSiC{i}=(-1/2)*inv(Ri)*Gsigmai*WaHSiC{i}
    muHSi=cell2mat(muHSiC);
    uSi = L_gi\(muHSi+FHi);
    muH_i = muHSiC{1};
    uSiC = mat2cell(uSi,mSi);
    u_i = uSiC{1};
%     VH_i = WcH_i'*sigmai;
    Q_i = e_i'*Qi*e_i;

    deltaH_i = WcH_i'*omega_i+Q_i+muH_i'*R_i*muH_i;
    pe_Wc_i=eta_c1i*Gamma_i*omega_i*deltaH_i/rho_i;
    pe_Wa_i=eta_c1i*GSigmai'*(R_i\GSigmai)*WaH_i*omega_i'*...
        WcH_i/(4*rho_i);
    M_i = M(i); % M is array containing number of extrapolation points 
    P_i = P(i);
    cl_Wc_i=0;
    cl_Wa_i=0;
    eSi_histC=e_histC{:,i};
    xSi_histC=x_histC{:,i};
    x0_k=x0_histC{i};
    parfor k=1:M_i
        eSikC=mat2cell(eSi_histC{k},n*ones(length(Si),1),1);
        xSikC=mat2cell(xSi_histC{k},n*ones(length(Si),1),1);
        e_i_k=eSikC{1};
%         x0_k= xSikC{1}-eSikC{1}; %Fix this!!!
        muHSikC = ADPCLNNTSmuSi(Si,xSikC,eSikC,mSi,NmSiC,ASi,A0Si,xDijSiC,RSiC,WaHSiC,Nset,n);
        omega_i_k = ADPCLNNTSOmegai(Si,xSikC,eSikC,thetaHSiC,WaHSiC,x0_k,mSi,NmSiC,ASi,A0Si,xDijSiC,xDi0SiC,RSiC,Nset,n,t);
        rho_i_k = 1+nu_i*omega_i_k'*Gamma_i*omega_i_k;
        muH_i_k = muHSikC{1};
        Q_i_k = e_i_k'*Qi*e_i_k;
        deltaH_i_k = WcH_i'*omega_i_k+Q_i_k+muH_i_k'*R_i*muH_i_k;
        cl_Wc_i = cl_Wc_i + omega_i_k*deltaH_i_k/rho_i_k;
        
        GSigmai_k = ADPCLNNTSGSigmai(Si,xSikC,eSikC,mSi,NmSiC,ASi,A0Si,xDijSiC,Nset,n);
        cl_Wa_i = cl_Wa_i+GSigmai_k'*(R_i\GSigmai_k)*WaH_i*...
            omega_i_k'*WcH_i/rho_i_k;
    end
    cl_Wc_i=eta_c2i*Gamma_i*cl_Wc_i/M_i;
    cl_Wa_i=eta_c2i*cl_Wa_i/(4*M_i);
    
    WcH_DotC{i}=-pe_Wc_i-cl_Wc_i;
    WaH_DotC{i}=-eta_a1i*(WaH_i-WcH_i)-eta_a2i*WaH_i+...
        pe_Wa_i+cl_Wa_i;
    Gamma_i_Dot=(beta_i*Gamma_i-eta_c1i*Gamma_i*(omega_i*omega_i'/rho_i^2)*...
        Gamma_i)*(norm(Gamma_i)<GammaBar_i);
    Gamma_DotC{i}=reshape(Gamma_i_Dot,(L(i))^2,1);
    
    % System ID
    k_i=kmat(i);
    GammaTheta_i=GammaThetaC{i};
    kTheta_i=kTheta(i);
    M_theta_i=M_theta(i);
%     u_i=(u_i+2*sin(2*t)*(t<=3));
    [xH_Dot_i,thetaH_Dot_i,SIDConLear_i,SIDRank_i]=...
        ADPCLNNTSSysID(i,xH_i,thetaH_i,u_i,t,...
        x_i,k_i,GammaTheta_i,kTheta_i,M_theta_i,fs,F,g);
    xH_DotC{i}=xH_Dot_i;
    vecthetaH_DotC{i}=reshape(thetaH_Dot_i,n*P(i),1);
    SIDConLearC{i}=SIDConLear_i;
    SIDRankC{i}=SIDRank_i;
    f_i = ADPCLNNTSf(i,x_i);
    g_i = ADPCLNNTSg(i,x_i);
%     pe=2*sin(exp(1)*t)+cos(pi*t);
    x_DotC{i} = f_i + g_i*(u_i);
    muC{i}=muH_i;
    uC{i}=u_i;
end
Zdot = [cell2mat(x_DotC);cell2mat(xH_DotC);...
    cell2mat(vecthetaH_DotC);...
    cell2mat(WcH_DotC);cell2mat(WaH_DotC);...
    cell2mat(Gamma_DotC)];
EC{1,round(t*1000+1)}=cell2mat(eC);
UC{1,round(t*1000+1)}=cell2mat(uC);
MUC{1,round(t*1000+1)}=cell2mat(muC);
fprintf('\b\b\b\b\b\b\b\b\b');   
end


