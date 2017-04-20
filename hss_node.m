% Author: Milinda Fernando
% School of Computing, University of Utah
% This file contains the fields that needed to perform HSS factorization of
% a matrix. 


classdef hss_node

    
    properties
        m_uiI;                  % Row ids of the partition
        m_uiJ;                  % Column ids of the partition. 
       
        m_uiRow;                % Row ids for the ID. 
        m_uiCol;                % Col ids for the ID
        
        m_uiB1;                 % size (k x k) core matrix child 1
        m_uiB2;                 % size (k x k) core matrix child 2
        m_uiD;                  % diagonal core. Only available to leaf nodes. 
        
        m_uiU;                  % ID matrix for the row sapce. 
        m_uiV;                  % ID matrix for the column sapce. 
        
        m_uiR_Row;              % Row part of the random matrix
        m_uiR_Col;              % Col part of the random matrix
        
        m_uiS_Row;              % Row part of the Sample matrix
        m_uiS_Col;              % Column part of the sample matrix
        
        
        m_uiY;                  % store the answer for matvec operations. 
        m_uib;                  % helper variable to compute the matvec. 
        
         
    end
    
    methods
        
        % default constructor. 
        function [obj]=hss_node(rowIds,colIds)
            obj.m_uiI=rowIds;
            obj.m_uiJ=colIds;
        end 
        
    end
    
end
