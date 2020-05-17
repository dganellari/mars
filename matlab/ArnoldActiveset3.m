


function [u,lambda,WorkingSet] = ArnoldActiveset3(A,B,b,c)

% we solve for min H, with H=0.5 x' A x - x' f - lambda (c-B x)
% structure of the problem
% |A B'| |x     |= |b|
% |B 0 | |lambda|  |c|
% in particular ST= matrix for quality constraints, CT= corresponding rhs
% in particular BT= matrix for inequality constraints, C= corresponding rhs

toll2=10^(-12);
toll=10^(-9);
n=length(b);
WorkingSet=zeros(n,1);
nconstraint=length(c);


continueplease=true;


%WorkingSet=[];
constraint_tot=[];

cont=0;
 while(continueplease) 
     
     w=find(WorkingSet>0);
     

B_k = B(w,:);
c_k = c(w);

LW=length(w);

M_k=[A    B_k';
     B_k  sparse(LW,LW);];


F_k=[b;c_k];


% use iterative refinement to compute the solution of the system:
% M_k * correction= F_k
% by doing 
% 1) M_k * xx= F_k
% 2) M_k * cc= F_k - M_k * xx
% 3) correction=xx+cc

[L,U,P] = lu(full(M_k));
yy=L\(P*F_k);
correction=U\yy;

for hh=1:1
res=F_k-M_k*correction;
yy=L\(P*res);
cc=U\yy;
correction=correction+cc;
end

%correction=M_k\F_k;

u_k=correction(1:n);

lambda_k=correction(1+n:end);
WorkingSetOld=WorkingSet;

res=b- A* u_k;

continueplease=false;
for i=1:n
    
    if(WorkingSet(i)>0 && res(i)<0)
        WorkingSet(i)=0;
        continueplease=true;      
    end
     
end

for i=1:n
    
    if(u_k(i)>c(i))
        WorkingSet(i)=1;
        continueplease=true;      
    end
     
end



end

 u=u_k;
 lambda=lambda_k;
 end





