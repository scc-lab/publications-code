function ei=ADPCLNNTSe(Si,xSiC,ASi,A0Si,xDSiC,xDi0,Nmi,x0)
xi = xSiC{1};
ei = A0Si(1,1)*(xi-x0-xDi0);
if isempty(Nmi)~=1
    for jt=1:length(Nmi)
        Node_jt_no = Nmi(jt); % agent number at jt'th position in Nmk
        j=Si==Node_jt_no; % logical array corresponding to position of 
                             % agent Node_jt_no in Si
        xj=xSiC{j}; % State of agent Node_jt_no
        xdij=xDSiC{1,j}; % Relative state of agent Node_jt_no
        ei=ei+ASi(1,j)*(xi-xj-xdij);
    end
end
