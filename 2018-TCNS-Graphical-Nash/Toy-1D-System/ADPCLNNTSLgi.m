function L_gi = ADPCLNNTSLgi(Si,xSiC,xDijSiC,ASi,A0Si,mSi,NmSiC)
L_giC=cell(length(Si),length(Si));
	for j=1:length(Si)
		for k=1:length(Si)
			if k==j
				L_giC{j,k}=(sum(ASi(j,:))+A0Si(j,j))*eye(mSi(j)); % note that sum over a row of ASi is the same as sum over the neighborhood set since all the other elements are zero.
			else
				xk=xSiC{k};
				xdjk=xDijSiC{j,k};
				node_j_no=Si(j);
				node_k_no=Si(k);
				gjk=ADPCLNNTSgplus(node_j_no,xk+xdjk)*ADPCLNNTSg(node_k_no,xk);
				L_giC{j,k}=-ASi(j,k)*gjk;
			end
		end
	end

	L_gi=cell2mat(L_giC);