function x=uzawa_patch_smoother2(Awithbc,bwithbc,x,A,b,l,u,constrained_dofs,d,RT2dofs,Multdofs,All2dofs,M,node_multiplier_dofs,Dirichlet,maxiter)


% C submatrix of speyee
% d contains bc

n=length(RT2dofs);

nA=length(u);

D=speye(nA);
cont=1;







%Dirichlet=[1,0,1,0,0,0,0,1,0,0,0,0,0]; 





tmp=0;

energy(cont)=0.5*x'*A*x-b'*x;
for it=1:maxiter
    it
    for i= 1:n     
        
        
        rt_dofs=RT2dofs{i};
        all_dofs=All2dofs{i};
        mult_dofs=Multdofs{i};
        
%         Ads=1:nA;
%         Bds=1+nA:length(b);
        
        
        
        
        n_rt_loc=length(rt_dofs);
        n_mult_loc=length(mult_dofs);
        n_all_loc=length(all_dofs);
        
        n_mult_loc=0;
        
        
        Aloc=full(Awithbc(all_dofs,all_dofs));

        
        bloc=bwithbc(all_dofs)-Awithbc(all_dofs,:)*x;
        
        
        
%         Mloc=full(M(mult_dofs,mult_dofs));
%         
%         Zeros=sparse(n_mult_loc,n_rt_loc);
%         ZerosMult=sparse(n_mult_loc,n_mult_loc);
        
        
      
        
        

        
        c_l=l(rt_dofs)-x(rt_dofs);
        c_u=u(rt_dofs)-x(rt_dofs);
        
        cloc=inf*ones(n_all_loc+n_mult_loc,1);
        cloc(1:n_rt_loc)=c_u;
        
        
%         Dloc=[Zeros Mloc];
        
      
        
%           Aloc=[Hloc Dloc'; Dloc ZerosMult];
        
%         bloc(end+1:end+n_mult_loc)=0;
                    
                
      
        
%         Atmp=[Aloc Bloc_new' ;Bloc_new zeros(r,r);]; 
%         btmp=[res_u;res_g_new;];
%         ctmp=c_u;
%         ctmp(1+nAloc:nAloc+r)=inf;
%        
%          
%          res_g=[res_g;zeros(n_bc,1)];
%          size_B=size(Bloc);
%          
%          tmp=min(tmp,rank(full(Bloc))-size_B(1));
         
%           Atmp=[Aloc Bloc' ;Bloc zeros(size_B(1),size_B(1));]; 
%           Atmp_size=size(Atmp);
         
% [i-1, rank(Aloc),size(Aloc,1)]
         
         
%          v(i,:)=[i,rank(full(Bloc))-size_B(1),rank(full(Atmp))-Atmp_size(1) ];
         
%          c = quadprog(Aloc,-res_u,[],[],Bloc,res_g_new,c_l,c_u);

         
         
% in this case, the problem is well posed
if( Dirichlet(i))
    
 [c,lambda,WorkingSet] = ArnoldActiveset3(Aloc,speye(n_all_loc+n_mult_loc),bloc,cloc);
 
% in this case we have to remove 3 rigid body motions or 1 displacament
% (if there is at least one condition on the stresses instead of two)
else
    
    only_disp=false;
    
     Dloc=full(Awithbc(rt_dofs,rt_dofs));
     
     for j=1:2:n_rt_loc
         % we consider condition on the first stress, so displacement can
         % be free as the first component
         if(Dloc(j,j)==1 && Dloc(j+1,j+1)~=1 && sum(Dloc(j,:))==1 && sum(Dloc(:,j))==1 && sum(Dloc(j+1,:))~=1 && sum(Dloc(:,j+1))~=1 )
             
             
             only_disp=true;
                                     
         end
         
         
     end
     
     
     if(0)%only_disp)
         
         all_dofs(n_rt_loc+1)=[];
         n_all_loc=length(all_dofs);
         Aloc(n_rt_loc+1,:)=[];
         Aloc(:,n_rt_loc+1)=[];
         bloc(1+n_rt_loc)=[];
         cloc(1+n_rt_loc)=[];       
         [c,lambda,WorkingSet] = ArnoldActiveset3(Aloc,speye(n_all_loc),bloc,cloc);
         
     else
         
         node_mult=node_multiplier_dofs{i};
         node_mult=node_mult([1,2,end]);
         
         for s=1:length(node_mult)
             
             all_dofs(find(node_mult(s)==all_dofs))=[];                    
         end
         n_all_loc=length(all_dofs);
         Aloc=full(Awithbc(all_dofs,all_dofs));  
         bloc=bwithbc(all_dofs)-Awithbc(all_dofs,:)*x;
         cloc=inf*ones(n_all_loc+n_mult_loc,1);
         cloc(1:n_rt_loc)=c_u;       
         [c,lambda,WorkingSet] = ArnoldActiveset3(Aloc,speye(n_all_loc),bloc,cloc);
         
     end
    
    
end

[it,i,norm(c)]
        x(all_dofs)=x(all_dofs)+c(1:n_all_loc);
        
        
%                  x(Ads)=x(Ads)+vvv(1:nAloc);
%         x(Bds(indipendent_rows))=x(Bds(indipendent_rows))+c(1+nAloc:end);
        
        cont=cont+1;
        

        energy(cont)=0.5*x'*A*x-b'*x;
        if(energy(cont)>energy(cont-1)+0.000000001)
            fermami=1;
            
        end
        
    end
    
end

tmp

figure
plot(energy);




end