function calGi = ADPCLNNTScalGi(xSiC,mSi,Si,NmSiC,ASi,A0Si,xDSiC)
	
	i=Si(1);
    xi=xSiC{1};
	L_gi = ADPCLNNTSLgi(Si,xSiC,xDSiC,ASi,A0Si,mSi);
	L_giC = mat2cell(inv(L_gi),mSi);
	calGi = ADPCLNNTSg(i,xi)*L_giC{1};
