function [xSi,x0] = ADPCLNNTSGenStateFromError(n,ASi,A0Si,ESi,xDijSiC,xDi0SiC)
DSi = diag(sum(ASi,2));
LSi = DSi-ASi;
si=size(ASi,1);
eSi=ESi(1:(end-n),1);
xi=ESi((end-n+1):end,1);

% Find xdi0 corresponding to the first agent connected to the leader
% Find index of the first agent connected to the leader
% xdi0Si=cell2mat(xDi0SiC);
k=find(diag(A0Si), 1, 'first'); % index of the first nonzero element
% k=ceil(k1/n); % index of the corresponding agent
xdk0 = xDi0SiC{k};

% Let xd10 = 1 and find all other xdi0
xDi0SiCrel=cell(si,1);
xDi0SiCrel{1}=zeros(size(xi));
for j=2:si
    xDi0SiCrel{j}=xDi0SiCrel{j-1}-xDijSiC{j-1,j};
end

% Find the difference between xdSi0(k) and xdk0 and add it to all the xdi0
d=xdk0-xDi0SiCrel{k};
xdi0Si=cell2mat(xDi0SiCrel)+kron(ones(si,1),d);

% ySi = xSi - xdi0Si - x0Si = kron((LSi+A0Si),eye(n))\eSi
ySi = kron((LSi+A0Si),eye(n))\eSi;

% Compute x0 = xi-xdi0-yi and x0Si
xdi0 = xdi0Si(1:n,1);
yi = ySi(1:n,1);
x0 = xi-xdi0-yi;
x0Si = kron(ones(si,1),x0);

% Compute xSi = ySi + xdi0Si + x0Si
xSi = ySi + xdi0Si + x0Si;
