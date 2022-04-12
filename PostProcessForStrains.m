function [epsilon_ave]=PostProcessForStrains(iele,nnpe,nqptv,wtqp,np,x,y,z,sfac,dndxi,dndet,dndze,u)


% Calculate strain tensors from the displacement field 
%  also computes average by integrating over the element

detj = zeros(1,nqptv);
dndx = zeros(nnpe,nqptv);
dndy = zeros(nnpe,nqptv);
dndz = zeros(nnpe,nqptv);

epsilon = zeros(6,nqptv);
epsilon_ave = zeros(6,1);
elevol = 0.0;
u_el = zeros(3*nnpe,1);

for j=1:1:nnpe
     
  j1=np(iele,j) ;  

  u_el(3*j-2)=u(3*j1-2);
  u_el(3*j-1)=u(3*j1-1);
  u_el(3*j)=u(3*j1);
end


[dndx,dndy,dndz,detj,elevol] = SFDerivatives(iele,nnpe,nqptv,wtqp,np,x,y,z,dndxi,dndet,dndze);

bigB = BigBmat(nnpe,nqptv,dndx,dndy,dndz);

% Strain from small strain matrix
  
    
for i=1:1:nqptv  
epsilon(:,i) =bigB(:,:,i)*u_el(:);
epsilon_ave(:) = epsilon_ave(:) + wtqp(i)*detj(i)*epsilon(:,i);
end

epsilon_ave(:) = epsilon_ave(:)/elevol;

if(epsilon_ave(1)>-epsilon_ave(2)/2)
  a=1;  
end



