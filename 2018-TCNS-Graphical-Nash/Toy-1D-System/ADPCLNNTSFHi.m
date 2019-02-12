function FHi=ADPCLNNTSFHi(xSiC,thetaHSiC,x0,Si,NmSiC,ASi,A0Si,xDSiC,xD0SiC,t)
	%Si: subgraph set with agent i in first position
	%NmSiC: Cell of neighborhood sets corresponding to agents in Si
	%mSi: Vector of dimensions of inputs corresponding to agents in Si
	%XDSiC: Cell of vectors xdij corresponding to agents in Si
	%XD0SiC: Cell of vectors xdi0 corresponding to agents in Si
	FHiC = cell(size(Si));
	for k=1:length(Si)
		Node_k_no = Si(k); % agent number at k'th position in Si 
		Nmk=NmSiC{k}; % Neighborhood set for agent at k'th position in Si
		xk=xSiC{k}; % State of the agent at k'th position in Si
		thetaHk=thetaHSiC{k}; % Weight estimate for the agent
											   % at k'th position in Si
		row=A0Si(k,k)*ADPCLNNTSfHij(Node_k_no,0,xk,thetaHk,...
			x0,1,xD0SiC{k},t);
		if isempty(Nmk)~=1
			for jt=1:length(Nmk)
				Node_jt_no = Nmk(jt); % agent number at jt'th position in Nmk
				j=Si==Node_jt_no; % logical array corresponding to position of 
								 % agent Node_l_no in Si
				xj=xSiC{j}; % State of agent Node_l_no
				thetaHj=thetaHSiC{j}; % Weight of agent Node_l_no
				row = row + ASi(k,j)*ADPCLNNTSfHij(Node_k_no,...
					Node_jt_no,xk,thetaHk,xj,thetaHj,xDSiC{k,j},t);
			end
		end
		FHiC{k,1} = row;
	end
	FHi=cell2mat(FHiC);