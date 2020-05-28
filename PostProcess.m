function [epsilon,epsilon_ave,sigma,sigma_ave,dsigma]=PostProcess(iele,rotated_rmatrix,nnpe,nqptv,wtqp,np,x,y,z,sfac,dndxi,dndet,dndze,u)


% Calculate the stress and strain tensors from the displacement field 
%  also computes average by integrating over the element

detj = zeros(1,nqptv);
dndx = zeros(nnpe,nqptv);
dndy = zeros(nnpe,nqptv);
dndz = zeros(nnpe,nqptv);

sigma = zeros(6,nqptv);
epsilon = zeros(6,nqptv);
dsigma = zeros(6,nqptv);
sigma_ave = zeros(6,1);
epsilon_ave = zeros(6,1);
elevol = 0.0;

for j=1:1:nnpe
     
  j1=np(iele,j) ;  

  u_el(3*j-2,1)=u(3*j1-2);
  u_el(3*j-1,1)=u(3*j1-1);
  u_el(3*j,1)=u(3*j1);
end


[dndx,dndy,dndz,detj,elevol] = SFDerivatives(iele,nnpe,nqptv,wtqp,np,x,y,z,dndxi,dndet,dndze);

bigB = BigBmat(nnpe,nqptv,dndx,dndy,dndz);

% Strain from small strain matrix
  
    
for i=1:1:nqptv  
    
epsilon_qp(:) =bigB(:,:,i)*u_el(:);

% Stress using linear elasticity
sigma_qp=rotated_rmatrix*epsilon_qp';

%deviatoric stress
trsigma = (sigma_qp(1)+sigma_qp(2)+sigma_qp(3))/3.;
meanstress = trsigma*[1;1;1;0;0;0];
dsigma_qp = sigma_qp-meanstress;

sigma(:,i) = sigma_qp(:);
dsigma(:,i) = dsigma_qp(:);
epsilon(:,i) = epsilon_qp(:);

sigma_ave(:) = sigma_ave(:)+ sigma_qp(:)*wtqp(i)*detj(i);
epsilon_ave(:) = epsilon_ave(:)+ epsilon_qp(:)*wtqp(i)*detj(i);

end

sigma_ave(:) = sigma_ave(:)/elevol;
epsilon_ave(:) = epsilon_ave(:)/elevol;



