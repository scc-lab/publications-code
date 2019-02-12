function omegai = ADPCLNNTSOmegai(Si,xSiC,eSiC,...
    thetaHSiC,WaHSiC,x0,mSi,NmSiC,ASi,...
    A0Si,xDijSiC,xDi0SiC,RSiC,NSet,n,t)

	omegai = 0;
	[SiC,ASiC,A0SiC,SiRelC,~]=ADPCLNNTSGenerateSubgraphs(ASi,A0Si,Si);
	i = Si(1);
	 for j=1:length(Si)
		 
		Node_j_no = Si(j); % agent number at j'th position in Si 
		
		Sj = SiC{j};
		SjRel = SiRelC{j};
		xSjC = xSiC(SjRel);
		eSjC = eSiC(SjRel);
		RSjC = RSiC(SjRel);
		thetaHSjC=thetaHSiC(SjRel);
		mSj=mSi(SjRel);
		NmSjC=NmSiC(SjRel);
		ASj=ASiC{j};
		A0Sj=A0SiC{j};
		xDjkSjC=xDijSiC(SjRel,SjRel);
		xDj0SjC=xDi0SiC(SjRel);
		WaHSjC = WaHSiC(SjRel);
		
		scrFHj = ADPCLNNTSscrFHi(xSjC,thetaHSjC,x0,mSj,Sj,NmSjC,ASj,A0Sj,xDjkSjC,xDj0SjC,t);
		calFHj = ADPCLNNTScalFHi(xSjC,thetaHSjC,x0,mSj,Sj,NmSjC,ASj,A0Sj,xDjkSjC,xDj0SjC,t);
		scrGj = ADPCLNNTSscrGi(xSjC,mSj,Sj,NmSjC,ASj,A0Sj,xDjkSjC);
		calGj = ADPCLNNTScalGi(xSjC,mSj,Sj,NmSjC,ASj,A0Sj,xDjkSjC);
		
		muHSjC = ADPCLNNTSmuSi(Sj,xSjC,eSjC,mSj,NmSjC,ASj,A0Sj,xDjkSjC,RSjC,WaHSjC,NSet,n);
		muHSj = cell2mat(muHSjC);

		Grad_ejSig_i = ADPCLNNTSGradSigma(i,Node_j_no,eSiC,xSiC{1},Si,NSet,n);

		if j==1
			Grad_xiSig_i = ADPCLNNTSGradSigma(i,0,eSiC,xSiC{1},Si,NSet,n);
			omegai = omegai + Grad_xiSig_i*(calFHj+calGj*muHSj);
		end
		
		omegai = omegai + Grad_ejSig_i*(scrFHj+scrGj*muHSj);
	 end