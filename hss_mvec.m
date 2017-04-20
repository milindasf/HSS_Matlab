function [y] = hss_mvec(hss_tree,x)
% Author: Milinda Fernando
% School of Computing , University of Utah
% email address: milinda@cs.utah.edu
% 15/04/2017

%hss_mvec - perform matrix vector multiplication using,  hss factorization
%  
% Assumptions: HSS factorization has already be performed. 
%
% Syntax:  [y] = hss_mvec(hss_tree,x)
%
% Inputs:
%    hss_tree - HSS factorization of the input matrix A
%    x - input vector to multiply with
%   
%
% Outputs:
%   y - output vector such that y=Ax
%
%
% Other m-files required: hss_node.m
% MAT-files required: node
%
% See also: 
%------------- BEGIN CODE --------------


numTotalNodes=size(hss_tree,2)+1;
L=int32(log2(numTotalNodes));
numLeafNodes=int32(2^(L-1));
N=size(x,1);

y=zeros(N,1);

for n=0:(numLeafNodes-1)
    parent=2^(L-1) +n;
    I_l=(idivide(n*N,numLeafNodes)+1):idivide((n+1)*N,numLeafNodes);
    hss_tree(parent).m_uiY = hss_tree(parent).m_uiV'*x(I_l);
end


for lev=(L-1):-1:1
    
    numNodes=int32(2^(lev-1));

    for n=0:(numNodes-1)
    
      parent=2^(lev-1) +n;
      child1=2^(lev) + 2*n;
      child2=2^(lev) + 2*n+1;  
      hss_tree(parent).m_uiY=hss_tree(parent).m_uiV'*vertcat(hss_tree(child1).m_uiY,hss_tree(child2).m_uiY);   
    
    end
    
end

rank_k=size(hss_tree(1).m_uiU,2);
hss_tree(1).m_uib=zeros(rank_k,1);

for lev=1:(L-1)
     numNodes=int32(2^(lev-1));

    for n=0:(numNodes-1)
      parent=2^(lev-1) +n;
      child1=2^(lev) + 2*n;
      child2=2^(lev) + 2*n+1;  
      
      W=hss_tree(parent).m_uiU*hss_tree(parent).m_uib;
     
      hss_tree(child1).m_uib = hss_tree(parent).m_uiB1*hss_tree(child2).m_uiY+ W(1:rank_k);
      hss_tree(child2).m_uib = hss_tree(parent).m_uiB2*hss_tree(child1).m_uiY+ W(rank_k+1:end);
      
      
    end
end



for n=0:(numLeafNodes-1)
    parent=2^(L-1) +n;
    I_l=(idivide(n*N,numLeafNodes)+1):idivide((n+1)*N,numLeafNodes);
    y(I_l) = hss_tree(parent).m_uiU*hss_tree(parent).m_uib + hss_tree(parent).m_uiD*x(I_l);
end



end