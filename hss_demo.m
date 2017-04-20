function demo(matSize, rank) 
% DEMO Access examples via Help browser. 
%
%   DEMO opens the Help browser to MATLAB Examples.
%
%   DEMO TYPE NAME opens the examples for the product matching NAME of
%   and TYPE, as defined in that product's info.xml or demos.xml
%   file.
%   
%   Examples:
%       demo 'matlab'
%       demo 'toolbox' 'signal'
%
%   See also DOC.


A=magic(matSize);
hss_tree=hss(A,rank,16);

x=rand(matSize,1);

tic;
y=A*x;
dense_mvec_time=toc;

tic;
y1=hss_mvec(hss_tree,x);
hss_mvec_time=toc;

fprintf('matlab matvec (s): %f \n',dense_mvec_time);
fprintf('hss mvec (s): %f \n ',hss_mvec_time);
fprintf('speed up : %f \n',(dense_mvec_time/hss_mvec_time))

relError=norm(y1-y)/norm(y);
fprintf('Rel Error : %f\n',relError);

end
