function muHSiC = ADPCLNNTSmuSi(Si,xSiC,eSiC,mSi,NmSiC,ASi,A0Si,xDijSiC,RSiC,WaHSiC,NSet,n)
	[SiC,ASiC,A0SiC,SiRelC,~]=ADPCLNNTSGenerateSubgraphs(ASi,A0Si,Si);
	muHSiC = cell(size(Si));
	for j=1:length(Si)
		
	%     Node_j_no = Si(j); % agent number at j'th position in Si
		Sj = SiC{j};
		SjRel = SiRelC{j};
		xSjC = xSiC(SjRel);
		mSj=mSi(SjRel);
		NmSjC=NmSiC(SjRel);
		ASj=ASiC{j};
		A0Sj=A0SiC{j};
		xDjkSjC=xDijSiC(SjRel,SjRel);
		Rj=RSiC{j};
		eSjC = eSiC(SjRel);    
		
		Gsigmaj = ADPCLNNTSGSigmai(Sj,xSjC,eSjC,mSj,NmSjC,ASj,A0Sj,xDjkSjC,NSet,n);

		muHSiC{j} = (-1/2)*inv(Rj)*Gsigmaj*WaHSiC{j};
	end