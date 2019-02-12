function GSigmai = ADPCLNNTSGSigmai(Si,xSiC,eSiC,mSi,NmSiC,ASi,A0Si,xDijSiC,NSet,n)
    
i=Si(1);
[SiC,ASiC,A0SiC,SiRelC,~]=ADPCLNNTSGenerateSubgraphs(ASi,A0Si,Si);
GSigmai = 0;

 for j=1:length(Si)

    Node_j_no = Si(j); % agent number at j'th position in Si 
    Grad_ejSig_i = ADPCLNNTSGradSigma(i,Node_j_no,eSiC,xSiC{1},Si,NSet,n);

    Sj = SiC{j};
    SjRel = SiRelC{j};
    xSjC = xSiC(SjRel);
    mSj=mSi(SjRel);
    NmSjC=NmSiC(SjRel);
    ASj=ASiC{j};
    A0Sj=A0SiC{j};
    xDjkSjC=xDijSiC(SjRel,SjRel);
    
    partialmuSjmui = zeros(sum(mSj),mSi(1));
    partialmuSjmuiC = mat2cell(partialmuSjmui,mSj);
    if nnz(Sj==i)~=0
        partialmuSjmuiC{Sj==i} = eye(mSi(1));
    end
    partialmuSjmui = cell2mat(partialmuSjmuiC);

    scrGj = ADPCLNNTSscrGi(xSjC,mSj,Sj,NmSjC,ASj,A0Sj,xDjkSjC);
	calGj = ADPCLNNTScalGi(xSjC,mSj,Sj,NmSjC,ASj,A0Sj,xDjkSjC);

    if j==1
        Grad_xiSig_i = ADPCLNNTSGradSigma(i,0,eSiC,xSiC{1},Si,NSet,n);
        GSigmai = GSigmai + partialmuSjmui'*calGj'*Grad_xiSig_i';
    end

    GSigmai = GSigmai + partialmuSjmui'*scrGj'*Grad_ejSig_i';
 end