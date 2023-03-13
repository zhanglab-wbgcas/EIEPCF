%LP_TGN_infer: Linear Programming Method for infering the TF-gene
%network,Y is the gene expression profile ,which is approximate to TF
%activity.X is the TF activity or gene(coding TF)expression profile
% degree_connect is the ratio of nonzeros stem  or the sparseness.
% parameter is the value below which will be deleted in the infered
% network,or called the error
function [J]=LP_TGN(Y,X,lamda)
[n,m]=size(Y);
[p,~]=size(X);
% lamda=1;
c=n*m; h=n*p;
f=[ones(1,2*c),lamda*ones(1,2*h)]';
Y_1=Y';
beq=Y_1(:);
A1=sparse(1:c,1:c,ones(1,c),c,c);
A2=sparse(1:c,1:c,-ones(1,c),c,c);
Z=X';
Z_1=Z;
for i=1:n-1
    Z_1=blkdiag(Z_1,Z);
    Z_1=sparse(Z_1);
   % i
end
Z_2=-Z_1;
Aeq=[A1,A2,Z_1,Z_2];
Aeq=sparse(Aeq);
clear A1 A2 Z_1 Z_2 Z Y_1;
lb=zeros(2*(c+h),1);
x=linprog(f,[],[],Aeq,beq,lb);
x=x';
s=zeros(n,p);
s=sparse(s);
t=zeros(n,p);
t=sparse(t);
J=zeros(n,p);
J=sparse(J);
for i=1:n
    for k=1:p
        s(i,k)=x(2*n*m+(i-1)*p+k);
        t(i,k)=x(2*n*m+n*p+(i-1)*p+k);
    end
end
% clear n p x Aeq beq lb f i k c h m q
J=s-t;
% clear s t  
J=full(J);
%J(abs(J)<parameter)=0;
end

%save J J
%dlmwrite('result.xls',J,'delimiter','\t','precision',8);
