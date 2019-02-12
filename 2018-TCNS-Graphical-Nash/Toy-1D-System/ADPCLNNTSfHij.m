function fHij = ADPCLNNTSfHij(i,j,~,thetaHi,xj,thetaHj,xdij,t)

[sigthj,~]=ADPCLNNTSsigmatheta(j,xj,t);
[sigthdi,~]=ADPCLNNTSsigmatheta(i,xj+xdij,t);
fHij=ADPCLNNTSgplus(i,xj+xdij)*(thetaHj'*sigthj-thetaHi'*sigthdi);