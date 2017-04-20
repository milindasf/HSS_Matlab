function [hss_tree] = hss(A,k,p)
%hss - Perform HSS factorization of the input matrix A
%  
% Assumptions: Assume that matrix A is Heirarchically semi-seperable, and
% rank  of off diagonal blocks of A < k. 
%
%
% Syntax:  [B,U,V] = hss(A,20,10)
%
% Inputs:
%    A - Square matrix A that that we are interested in computing HSS
%     factors. 
%    k - Upper bound for the rank of off-diagonal blocks of k. 
%    p - oversampling parameter used for interpolative decomposition. 
%
% Outputs:
%    B - List of B matrices of all levels in the order of bottom up size(kxk). 
%    U - List of U matrices of all levels in the order of bottom up size (n/2^l xk)
%    V - List of V matrices of all levels in the order of bottom up. size (n/2^l xk)
% Example: 
%    Line 1 of example
%    Line 2 of example
%    Line 3 of example
%
% Other m-files required: hss_node.m
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Author: Milinda Fernando
% School of Computing , University of Utah
% email address: milinda@cs.utah.edu
% 15/04/2017

%------------- BEGIN CODE --------------

if(size(A,1)~=size(A,2))
    fprintf('Input matrix should be a square matrix.')
    return ;
end
N=int32(size(A,1));
R_Row=normrnd(0,1,N,(k+p)); % generates a Gaussian random matrix. 
R_Col=normrnd(0,1,N,(k+p)); % generates a Gaussian random matrix. 

SRow=A*R_Row; 
SCol=A*R_Col;

L=int32(log2(double((N/k))))+1;

% initialize HSS tree. 
nCount=1;
for lev=1:L
    for n=1:(2^(lev-1))
        hss_tree(nCount)=hss_node(1,1);
        nCount=nCount+1;
    end
end
%size(hss_tree)

for lev=L:-1:1
    
    %fprintf('lev: %i\n',lev);
    if(lev==L) % leaf level nodes. 
        numNodes=int32(2^(lev-1));
        for n=0:(numNodes-1)
            %fprintf('n: %d begin : %d end: %d \n',n,(idivide(n*N,numNodes)+1),idivide((n+1)*N,numNodes));
            I_l=(idivide(n*N,numNodes)+1):idivide((n+1)*N,numNodes);
            
            IRow_l=I_l;
            JCol_l=I_l;
            parent=2^(lev-1) +n;
            hss_tree(parent)=hss_node(IRow_l,JCol_l);
            hss_tree(parent).m_uiD=A(I_l,I_l);
            %disp(I_l);
            
            RRow_l=R_Row(I_l,:);
            RCol_l=R_Col(I_l,:);
            
            SRow_l=SRow(I_l,:)-A(I_l,I_l)*RRow_l;
            SCol_l=SCol(I_l,:)-A(I_l,I_l)*RCol_l;
            
            II=eye(k);
            [SKRow_l,RDRow_l,TRow_l]=id(SRow_l',k);
            [SKCol_l,RDCol_l,TCol_l]=id(SCol_l',k);
            
            hss_tree(parent).m_uiU(:,RDRow_l)=TRow_l;
            hss_tree(parent).m_uiU(:,SKRow_l)=II;
            hss_tree(parent).m_uiV(:,RDCol_l)=TCol_l;
            hss_tree(parent).m_uiV(:,SKCol_l)=II;
            
            hss_tree(parent).m_uiR_Row=hss_tree(parent).m_uiV'*RRow_l;
            hss_tree(parent).m_uiR_Col=hss_tree(parent).m_uiU'*RCol_l;
            hss_tree(parent).m_uiS_Row=SRow_l(SKRow_l,:);
            hss_tree(parent).m_uiS_Col=SCol_l(SKCol_l,:);
            hss_tree(parent).m_uiRow=IRow_l(SKRow_l);
            hss_tree(parent).m_uiCol=JCol_l(SKCol_l);
            
        end
        
    else
         numNodes=int32(2^(lev-1));
         for n=0:(numNodes-1)
            
            parent=2^(lev-1) +n;
            child1=2^(lev) + 2*n;
            child2=2^(lev) + 2*n+1;
            
            %fprintf('n: %d parent: %d, child1: %d child2:%d \n',n,parent,child1,child2);
            I_l=(idivide(n*N,numNodes)+1):idivide((n+1)*N,numNodes);
            
            I_v1_row=hss_tree(child1).m_uiRow;   I_v1_col=hss_tree(child1).m_uiCol;
            I_v2_row=hss_tree(child2).m_uiRow;   I_v2_col=hss_tree(child2).m_uiCol; 
            
            IRow_l=horzcat(I_v1_row,I_v2_row);
            JCol_l=horzcat(I_v1_col,I_v2_col);
            
            hss_tree(parent)=hss_node(I_l,I_l); 
            %disp(I_l);
        
            RRow_l=vertcat(hss_tree(child1).m_uiR_Row,hss_tree(child2).m_uiR_Row);
            RCol_l=vertcat(hss_tree(child1).m_uiR_Col,hss_tree(child2).m_uiR_Col);
            
            SRow_l=vertcat((hss_tree(child1).m_uiS_Row-A(I_v1_row,I_v2_col)*hss_tree(child2).m_uiR_Row),(hss_tree(child2).m_uiS_Row-A(I_v2_row,I_v1_col)*hss_tree(child1).m_uiR_Row));
            SCol_l=vertcat((hss_tree(child1).m_uiS_Col-A(I_v1_row,I_v2_col)*hss_tree(child2).m_uiR_Col),(hss_tree(child2).m_uiS_Col-A(I_v2_row,I_v1_col)*hss_tree(child1).m_uiR_Col));
                        
            II=eye(k);
            [SKRow_l,RDRow_l,TRow_l]=id(SRow_l',k);
            [SKCol_l,RDCol_l,TCol_l]=id(SCol_l',k);
            
            hss_tree(parent).m_uiU(:,RDRow_l)=TRow_l;
            hss_tree(parent).m_uiU(:,SKRow_l)=II;
            hss_tree(parent).m_uiV(:,RDCol_l)=TCol_l;
            hss_tree(parent).m_uiV(:,SKCol_l)=II;
            
            hss_tree(parent).m_uiV=hss_tree(parent).m_uiV';
            hss_tree(parent).m_uiU=hss_tree(parent).m_uiU';
            
            
            hss_tree(parent).m_uiR_Row=hss_tree(parent).m_uiV'*RRow_l;
            hss_tree(parent).m_uiR_Col=hss_tree(parent).m_uiU'*RCol_l;
            
            hss_tree(parent).m_uiS_Row=SRow_l(SKRow_l,:);
            hss_tree(parent).m_uiS_Col=SCol_l(SKCol_l,:);
            
            hss_tree(parent).m_uiRow=IRow_l(SKRow_l);
            hss_tree(parent).m_uiCol=JCol_l(SKCol_l);
            
            hss_tree(parent).m_uiB1=A(I_v1_row,I_v2_col);
            hss_tree(parent).m_uiB2=A(I_v2_row,I_v1_col);
            

         end
    end
    
    
    
end


%------------- END OF CODE -------------

end