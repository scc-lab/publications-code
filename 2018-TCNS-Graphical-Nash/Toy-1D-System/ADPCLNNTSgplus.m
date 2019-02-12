function gplusi = ADPCLNNTSgplus(i,xi)
if i==0
%    gplusi=[0 0];
	gplusi=0;
else
%     gplusi = [0 1/cos(2*(xi(1)))+2];
    gplusi=1/(cos(2*xi)+2);
end