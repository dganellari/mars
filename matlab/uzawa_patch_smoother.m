function x=uzawa_patch_smoother(Awithbc,bwithbc,x,A,b,l,u,constrained_dofs,d,A2dofs,B2dofs,maxiter)


% C submatrix of speyee
% d contains bc

n=length(A2dofs);

nA=length(u);

D=speye(nA);
cont=0;













tmp=0;
for it=1:maxiter
    
    for i=1:n     
        
        
        Ads=A2dofs{i};
        Bds=B2dofs{i};
        
%         Ads=1:nA;
%         Bds=1+nA:length(b);
        
        
        
        
        nAloc=length(Ads);
        nBloc=length(Bds);
        
        
        Aloc=full(Awithbc(Ads,Ads));
        Bloc=full(Awithbc(Bds,Ads));
        
        res_u=bwithbc(Ads)-Awithbc(Ads,:)*x;
        res_g=bwithbc(Bds)-Awithbc(Bds,:)*x;        
                
         
        
        r=1;
        are_indipendent_rows=zeros(nBloc,1);
        are_indipendent_rows(1)=1;
        for j=2:nBloc
            
            if(rank(Bloc(1:j,:))>r)
                r=r+1;
                are_indipendent_rows(j)=1;
            else
                are_indipendent_rows(j)=0;
            end                           
        end
        
        indipendent_rows=find(are_indipendent_rows>0);
        
        Bloc_new=Bloc;
        res_g_new=res_g;
        
        
        
        
%         Bloc_new=Bloc(indipendent_rows,:);
%         res_g_new=res_g(indipendent_rows);
        
        

        
        c_l=l(Ads)-x(Ads);
        c_u=u(Ads)-x(Ads);
        
      
        
        
        
        Atmp=[Aloc Bloc_new' ;Bloc_new zeros(r,r);]; 
        btmp=[res_u;res_g_new;];
        ctmp=c_u;
        ctmp(1+nAloc:nAloc+r)=inf;
        [c,lambda,WorkingSet] = ArnoldActiveset3(Atmp,speye(nAloc+r),btmp,ctmp);
%          
%          res_g=[res_g;zeros(n_bc,1)];
%          size_B=size(Bloc);
%          
%          tmp=min(tmp,rank(full(Bloc))-size_B(1));
         
%           Atmp=[Aloc Bloc' ;Bloc zeros(size_B(1),size_B(1));]; 
%           Atmp_size=size(Atmp);
         
         
         
%          v(i,:)=[i,rank(full(Bloc))-size_B(1),rank(full(Atmp))-Atmp_size(1) ];
         
%          c = quadprog(Aloc,-res_u,[],[],Bloc,res_g_new,c_l,c_u);
if(norm(c)>0.1)
    
    
    
    fermami=1
    
end
      [i,(norm(btmp))]   
         
         
%          x(Ads)=x(Ads)+vvv(1:nAloc);
        x(Ads)=x(Ads)+c(1:nAloc);
        x(Bds)=x(Bds)+c(1+nAloc:end);
%         x(Bds(indipendent_rows))=x(Bds(indipendent_rows))+c(1+nAloc:end);
        
        cont=cont+1;
        
        energy(cont)=0.5*x'*A*x-b'*x;
        
    end
    
end

tmp

figure
plot(energy);




end