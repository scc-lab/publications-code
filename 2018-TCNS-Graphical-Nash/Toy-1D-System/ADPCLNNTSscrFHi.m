function scrFHi = ADPCLNNTSscrFHi(xSiC,thetaHSiC,x0,mSi,Si,NmSiC,ASi,A0Si,xDSiC,xD0SiC,t)

i=Si(1);
L_gi = ADPCLNNTSLgi(Si,xSiC,xDSiC,ASi,A0Si,mSi);
L_giC = mat2cell(inv(L_gi),mSi);
Nmi=NmSiC{1}; % Neighborhood set for agent i
xi=xSiC{1}; % State of agent i
thetaHi=thetaHSiC{1}; % Weight of agent i
[f0,~]=ADPCLNNTSsigmatheta(0,x0,t);
[sigthi,~]=ADPCLNNTSsigmatheta(i,xi,t);

FHi=ADPCLNNTSFHi(xSiC,thetaHSiC,x0,Si,NmSiC,ASi,A0Si,xDSiC,xD0SiC,t);

scrFHi = A0Si(1,1)*(thetaHi'*sigthi-f0)+A0Si(1,1)*ADPCLNNTSg(i,xi)*L_giC{1}*FHi;

if isempty(Nmi)~=1
    for jte=1:length(Nmi)
        Node_jte_no=Nmi(jte); % agent number at jte'th position in Nmi
        j=Si==Node_jte_no; % logical array corresponding to position of 
                         % agent Node_l_no in Si
        xj=xSiC{j}; % State of agent Node_l_no in Si
        thetaHj=thetaHSiC{j}; % Weight of agent Node_l_no in Si
        [sigthj,~]=ADPCLNNTSsigmatheta(Node_jte_no,xj,t);
        scrFHi = scrFHi + ASi(1,j)*...
            (thetaHi'*sigthi - thetaHj'*sigthj)...
            + ASi(1,j)*(ADPCLNNTSg(i,xi)*L_giC{1} - ...
            ADPCLNNTSg(Node_jte_no,xj)*L_giC{j})*FHi;
    end
end