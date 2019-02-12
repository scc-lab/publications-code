function scrGi=ADPCLNNTSscrGi(xSiC,mSi,Si,NmSiC,ASi,A0Si,xDSiC)
i=Si(1);
xi=xSiC{1};
L_gi = ADPCLNNTSLgi(Si,xSiC,xDSiC,ASi,A0Si,mSi);
L_giC = mat2cell(inv(L_gi),mSi);
scrGi = A0Si(1,1)*ADPCLNNTSg(i,xi)*L_giC{1};
Nmi=NmSiC{1};
if isempty(Nmi)~=1
    for jte=1:length(Nmi)
        Node_jte_no=Nmi(jte); % agent number at jte'th position in Nmi
        j=Si==Node_jte_no; % logical array corresponding to position of 
                         % agent Node_l_no in Si
        xj=xSiC{j}; % State of agent Node_l_no in Si
        scrGi = scrGi+ASi(1,j)*(ADPCLNNTSg(i,xi)*L_giC{1}-ADPCLNNTSg(Node_jte_no,xj)*L_giC{j});
    end
end