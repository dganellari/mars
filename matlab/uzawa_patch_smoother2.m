function [x,energy]=uzawa_patch_smoother2(energy,Awithbc,bwithbc,x,A,b,l,u,constrained_dofs,d,RT2dofs,Multdofs,All2dofs,RT2alldofs,Disp2alldofs,M,node_multiplier_dofs,Dirichlet,b_rbm,maxiter)


% C submatrix of speyee
% d contains bc

n=length(RT2dofs);

nA=length(u);

D=speye(nA);
cont=length(energy);







%Dirichlet=[1,0,1,0,0,0,0,1,0,0,0,0,0]; 





tmp=0;

energy(cont+1)=0.5*x'*A*x-b'*x;
for it=1:maxiter
    it
    cont_dirichlet=0;
    cont_contact=0;
    cont_other=0;
    for i= 1:n     
        
        
        rt_dofs=RT2dofs{i};
        rt_all_dofs=RT2alldofs{i};
        all_dofs=All2dofs{i};
        mult_dofs=Multdofs{i};
        disp_dofs=Disp2alldofs{i};
        
%         disp_dofs([end-1,end])=[];
        
        
%         tmp=setdiff(rt_all_dofs,rt_dofs);
%         rt_dofs=[rt_dofs;tmp];%([1,2,3,4,5,6,7,8])];
%         all_dofs=[all_dofs;tmp];%([1,2,3,4,5,6,7,8])];
        
        
%         Ads=1:nA;
%         Bds=1+nA:length(b);
        
        
        
        
        n_rt_loc=length(rt_dofs);
        n_mult_loc=length(mult_dofs);
        n_all_loc=length(all_dofs);
        n_disp_loc=length(disp_dofs);
        
        
%         n_mult_loc=0;
        
        
        Aloc=full(Awithbc(all_dofs,all_dofs));

        
        bloc=bwithbc(all_dofs)-Awithbc(all_dofs,:)*x;
        
        
        
%         Mloc=full(M(mult_dofs,mult_dofs));
%         
%         Zeros=sparse(n_mult_loc,n_rt_loc);
%         ZerosMult=sparse(n_mult_loc,n_mult_loc);
        
        
      
        
        

        
        c_l=l(rt_dofs)-x(rt_dofs);
        c_u=u(rt_dofs)-x(rt_dofs);
        
        cloc=inf*ones(n_all_loc,1);
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
    
 [c,lambda,WorkingSet] = ArnoldActiveset3(Aloc,speye(n_all_loc),bloc,cloc);
 
 cont_dirichlet=cont_dirichlet+1;
 
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
             break;
                                     
         end
         
         
     end
     
     
     if(0)%only_disp)
         cont_contact=cont_contact+1;
         
         
%          C=[];
%          
%          
%          C=b_rbm(disp_dofs)';
%          for j=2:2:length(C)
%              C(j)=0;           
%          end
%          
%          C=[zeros(1,n_rt_loc),C,zeros(1,n_all_loc-n_disp_loc-n_rt_loc)];
%          
%          Aloc=[Aloc C';
%                C 0];
%          
%          bloc(end+1)=0;
%          cloc(end+1)=inf;      
         
         b_aux=[0;0;0];
         c_aux=[inf;inf;inf];
         C=[b_rbm(mult_dofs)';
            b_rbm(mult_dofs)';
            b_rbm(mult_dofs)';];

         for j=2:2:length(disp_dofs)
             C(1,j)=0;       
             C(2,j-1)=0;       
         end
         for j=length(disp_dofs)+1:length(C)
             C(1,j)=0;
             C(2,j)=0;
         end
         
         for j=1:length(disp_dofs)
             C(3,j)=0;         
         end
         
         
         C=[zeros(3,n_rt_loc),C];
         
         
         [c,lambda,WorkingSet] = ArnoldActiveset4(Aloc,speye(n_all_loc),bloc,cloc,C);
         
     else
         cont_other=cont_other+1;
         
         C=[b_rbm(mult_dofs)';
            b_rbm(mult_dofs)';
            b_rbm(mult_dofs)';];

         for j=2:2:length(disp_dofs)
             C(1,j)=0;       
             C(2,j-1)=0;       
         end
         for j=length(disp_dofs)+1:length(C)
             C(1,j)=0;
             C(2,j)=0;
         end
         
         for j=1:length(disp_dofs)
             C(3,j)=0;         
         end
         
         
         C=[zeros(3,n_rt_loc),C];
         
           Aloc=[Aloc C';
               C zeros(3,3)];
         

         bloc(end+1:end+3)=0;
         cloc(end+1:end+3)=inf;       
         [c,lambda,WorkingSet] = ArnoldActiveset3(Aloc,speye(n_all_loc+3),bloc,cloc);
       

         
         
%          
%          
%          Mloc=full(M(disp_dofs,disp_dofs));
%          
%          Zeros=sparse(n_disp_loc,n_rt_loc);
%          
%          Zeros2=sparse(n_disp_loc,n_rt_loc);
%          
%          ZerosMult=sparse(n_disp_loc,n_disp_loc);
%          
%         Dloc=[Zeros Mloc ];
%         
%         
%          Dloc=[Dloc sparse(n_disp_loc, n_all_loc-n_disp_loc-n_rt_loc)];
%         
%         
%               
%         Aloc=[Aloc Dloc'; Dloc ZerosMult];
%         
%         bloc(end+1:end+n_disp_loc)=0;
%         
%         cloc(end+1:end+n_disp_loc)=inf;
        
        
        
         
         
%         Mloc=full(M(mult_dofs,mult_dofs));
%         
%         Zeros=sparse(n_mult_loc,n_rt_loc);
%         
%         ZerosMult=sparse(n_mult_loc,n_mult_loc);
% 
%         Dloc=[Zeros Mloc];
%               
%         Aloc=[Aloc Dloc'; Dloc ZerosMult];
%         
%         bloc(end+1:end+n_mult_loc)=0;
%         
%         cloc(end+1:end+n_mult_loc)=inf;
        
%         [c,lambda,WorkingSet] = ArnoldActiveset3(Aloc,speye(n_all_loc+n_disp_loc),bloc,cloc);
                    
                
      
   
%          node_mult=node_multiplier_dofs{i};
%          node_mult=node_mult([1,2,end]);
%          
%          tmp=setdiff(Multdofs{i},node_mult);
%          
%          node_mult=tmp([1,2,end]);
%          
%          for s=1:length(node_mult)
%              
%              all_dofs(find(node_mult(s)==all_dofs))=[];                    
%          end
%          n_all_loc=length(all_dofs);
%          Aloc=full(Awithbc(all_dofs,all_dofs));  
%          bloc=bwithbc(all_dofs)-Awithbc(all_dofs,:)*x;
%          cloc=inf*ones(n_all_loc+n_mult_loc,1);
%          cloc(1:n_rt_loc)=c_u;       
%          [c,lambda,WorkingSet] = ArnoldActiveset3(Aloc,speye(n_all_loc),bloc,cloc);
         
     end
    
    
end

[it,i,norm(c),min(c),    cont_dirichlet,   cont_contact, cont_other]
        x(all_dofs)=x(all_dofs)+c(1:n_all_loc);
        
%  dlmwrite('sol.txt',[length(x);x])
     
%                  x(Ads)=x(Ads)+vvv(1:nAloc);
%         x(Bds(indipendent_rows))=x(Bds(indipendent_rows))+c(1+nAloc:end);
        
        
        
        cont=cont+1;
        energy(cont)=0.5*x'*A*x-b'*x;
        

        
    end

    
end

tmp

figure
plot(energy);




end