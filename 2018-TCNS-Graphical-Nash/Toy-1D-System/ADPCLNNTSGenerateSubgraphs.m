% A = [0 0 1 0 0;0 0 0 0 0;1 0 0 0 0;1 0 0 0 0;0 1 0 0 0];
% A0 = [1 0 0 0 0;0 1 0 0 0;0 0 0 0 0;0 0 0 0 0;0 0 0 0 0];
% D = diag(sum(A,2));
% L = D-A;
% SiC=cell(size(A,1),1);
% ASC=cell(size(A,1),1);
% A0SC=cell(size(A,1),1);

function [SC,ASC,A0SC,SRelC,NmC] = ...
    ADPCLNNTSGenerateSubgraphs(A,A0,S)
    SC=cell(size(A,1),1);
    SRelC=cell(size(A,1),1);
    ASC=cell(size(A,1),1);
    A0SC=cell(size(A,1),1);
    NmC=cell(size(A,1),1);
    for i=1:size(A,1) 
    NmC{i}=S(A(i,:)~=0);% traverse rows index i belongs to the i'th row
    j=i;
    iter=1;
    Si=i; % set the first element of Si to i
    while(1)
    temp = find(A(j,:)~=0); % find neighbors of the j'th agent
%     NmC{i}=S(temp);
    if isempty(temp)
        % If no neighbors are found then we want to break the loop.
        % However, since all neighbors found in the first iteration are
        % added to Si, if we immediately break loop here, then we may 
        % miss extended neighbors of some agents in Si. So, we break only
        % if the number of while loop iterations is equal to the length of
        % Si. This means tH all the neighbots of agent i found in the
        % first iteration are checked for extended neighbors.
        if iter==length(Si) 
            break;
        else
            iter=iter+1;
            j=Si(iter);
            continue;
        end
    end
    ExtNeighFound=0;
    for k=1:length(temp)
        if nnz(Si==temp(k))==0
            Si=[Si; temp(k)]; %#ok<AGROW>
            ExtNeighFound=1;
        end
    end
    if ExtNeighFound==0;
        % same logic as before. If no extended neighbors ara found and if
        % the while loop iteration number equals length of Si (indicating
        % tH all the agents in Si have been checked for extended
        % neighbors), we exit the loop.
        if iter==length(Si)
            break;    
        else
            iter=iter+1;
            j=Si(iter);
            continue;
        end
    end
    iter=iter+1;
    j=Si(iter);
    end
    SRelC{i}=Si;
    SC{i}=S(Si);
    ASC{i}=A(Si,Si);
    A0SC{i}=A0(Si,Si);
    end
    

