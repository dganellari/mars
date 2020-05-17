clear all
clc
load ALnobc.dat
load ALwithbc.dat
load Abnobc.dat
load Abwithbc.dat
load Aconstrained_dofs.dat
load Aconstrained_vec.dat
load Aactive_x.dat
load Ac_sup.dat
load Ac_inf.dat
load AA2dofs.dat
load AB2dofs.dat
load AA2elements_dofs.dat
load AB2elements_dofs.dat

load Aall2dofs.dat
load ART2dofs.dat

load AAmultipliers.dat
load Amultipliersdofs.dat

load  Ais_patch_dirichlet.dat



n=length(Abnobc);
constrained_dofs=find(Aconstrained_dofs>0);
constrained_vec=Aconstrained_vec(constrained_dofs);
Lconstrained=length(constrained_dofs);


cont=0;
i=1;
while(i<length(AA2dofs))
    cont=cont+1;
    A2dofs{cont}=AA2dofs(i+1:i+AA2dofs(i));
    i=i+AA2dofs(i)+1;
    
end
cont=0;
i=1;
while(i<length(AB2elements_dofs))
    cont=cont+1;
    B2dofs{cont}=AB2elements_dofs(i+1:i+AB2elements_dofs(i));
    i=i+AB2elements_dofs(i)+1;
    
end

cont=0;
i=1;
while(i<length(AB2dofs))
    cont=cont+1;
    
    
    C2dofs{cont}=AB2dofs(i+1:i+AB2dofs(i));
    i=i+AB2dofs(i)+1;
    
end

cont=0;
i=1;
while(i<length(AA2elements_dofs))
    cont=cont+1;
    
    
    D2dofs{cont}=AA2elements_dofs(i+1:i+AA2elements_dofs(i));
    i=i+AA2elements_dofs(i)+1;
    
end










cont=0;
i=1;
while(i<length(Amultipliersdofs))
    cont=cont+1;
    
    
    Y2dofs{cont}=Amultipliersdofs(i+1:i+Amultipliersdofs(i));
    i=i+Amultipliersdofs(i)+1;
    
end


cont=0;
i=1;
while(i<length(Aall2dofs))
    cont=cont+1;
    
    
    Z2dofs{cont}=Aall2dofs(i+1:i+Aall2dofs(i));
    i=i+Aall2dofs(i)+1;
    
end

cont=0;
i=1;
while(i<length(ART2dofs))
    cont=cont+1;
    
    
    X2dofs{cont}=ART2dofs(i+1:i+ART2dofs(i));
    i=i+ART2dofs(i)+1;
    
end



for i=1:length(C2dofs)
    C2dofs{i}([1,2])=[];
    
    G2dofs{i}=setdiff(B2dofs{i},C2dofs{i});
    
end








cont=0;
i=1;
while(i<length(Ais_patch_dirichlet))
    cont=cont+1;
    
    
    is_patch_dirichlet{cont}=Ais_patch_dirichlet(i+1:i+Ais_patch_dirichlet(i));
    i=i+Ais_patch_dirichlet(i)+1;
    
end






A=spconvert(ALnobc);
Abc=spconvert(ALwithbc);
M=spconvert(AAmultipliers);
bbc=Abwithbc;
L=length(find(0~=diag(A)));
A1=A(1:L,1:L);



b1=Abnobc(1:L);
g1=Abnobc(1+L:end);

u1=Ac_sup(1:L);
l1=Ac_inf(1:L);



B1=sparse(n-L+Lconstrained,L);
D=sparse(Lconstrained,L);



B1(1:n-L,1:L)=A(L+1:end,1:L);
Htmp=sparse(Lconstrained,n);

for i=1:Lconstrained
    B1(n-L+i,constrained_dofs(i))=1;
    g1(n-L+i)=constrained_vec(i);
    Htmp(i,constrained_dofs(i))=1;
    
   
end


y = quadprog(A1,-b1,[],[],B1,g1,l1,u1);


Atmp=[A Htmp';Htmp sparse(Lconstrained,Lconstrained)];
btmp=[Abnobc;constrained_vec];
ctmp=Ac_sup;
ctmp(n+1:length(btmp))=inf;
%  [z1,z2,WorkingSet] = ArnoldActiveset3(Atmp,speye(length(btmp)),btmp,ctmp);

% y=Aactive_x(1:L);
% norm(x-y)
% 0.5*y'*A1*y-b1'*y
% 0.5*x'*A1*x-b1'*x


C1=sparse(Lconstrained,n);
C1(1:Lconstrained,1:L)=B1(n-L+1:end,1:L);
h1=g1(n-L+1:end);


maxiter=2000;
D=speye(length(b1));
d=zeros(L,1);
d(constrained_dofs)=constrained_vec;


b=Abnobc;
x=zeros(n,1);
 x=uzawa_patch_smoother2(Abc,bbc,x,A,b,l1,u1,constrained_dofs,d,X2dofs,Y2dofs,Z2dofs,M,C2dofs,is_patch_dirichlet{length(is_patch_dirichlet)},maxiter);
 
 
% y=uzawa_patch_smoother(Abc,bbc,zeros(n,1),A,b,l1,u1,constrained_dofs,d,D2dofs,B2dofs,maxiter)

dlmwrite('sol.txt',[length(y);y])
dlmwrite('sol.txt',[length(x);x])
% x = quadprog(A,-Abnobc,[],[],C1,h1,Ac_inf,Ac_sup);dlmwrite('sol.txt',[length(x);x])

