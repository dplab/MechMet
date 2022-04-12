function [epsilon_ave] = ComputeFELatticeStrains(numel,numnp,np,x,y,z,ngrains,grains,phases,shafac_structure,grain_rmatrix,bccode,bcvalue,Fsr,Fsv,fc,nomdisplacement)
%
% computes the grain-average lattice strain over finite element mesh 
%          

sfac = shafac_structure.sfac;
dndxi = shafac_structure.dndxi;
dndet = shafac_structure.dndet;
dndze = shafac_structure.dndze;
nqptv = shafac_structure.nqpt;
wtqp = shafac_structure.wtqp;
nnpe = shafac_structure.nnpe;

%
% allocate arrays for the global stiffness and force matrices
%  using sparse matrix solver 


Ksr=zeros(numel*9*nnpe*nnpe,1); 
Ksc=zeros(numel*9*nnpe*nnpe,1); 
Ksv=zeros(numel*9*nnpe*nnpe,1);
kc=1;

%
%defines body forces to be zero 
%(potentially could be used to include residual stresses)
%
gamvec = [ 0 0 0 ];
element_volumes = zeros(1,numel);
%
%
% Compute the elemental stiffness and force arrays and assemble
%
for   iele =1:1:numel
    
  
  [dndx,dndy,dndz,detj,elvolume] = SFDerivatives(iele,nnpe,nqptv,wtqp,np,x,y,z,dndxi,dndet,dndze);
  
  element_volumes(iele) = elvolume;
  
  bigB = BigBmat(nnpe,nqptv,dndx,dndy,dndz);
  
  bigN = BigNmat(nnpe,nqptv,sfac);
  
  igrain = grains(iele);
  rotated_rmatrix = grain_rmatrix(:,:,igrain); 
 
  se = ElementStiffnessMatrix(nnpe,nqptv,wtqp,bigB,detj,rotated_rmatrix);

  fe = ElementForceMatrix(nnpe,nqptv,wtqp,bigN,detj,gamvec);

  [Ksr,Ksc,Ksv,Fsr,Fsv,fc,kc] = AssembleSparse(iele,nnpe,np,se,fe,Ksr,Ksc,Ksv,Fsr,Fsv,fc,kc);

end

%
% Define sparse matrices and clean up
%
 f=sparse(Fsr,1,Fsv);
 sk=sparse(Ksr,Ksc,Ksv);
 clear Frr Fsv Ksr Ksc Ksv;
%   
% 
%  Solve for displacements using sparse direct solver 
%    and return result in a full vector
%
% Apply boundary conditions
% 
  [sk,f] = BCSparse(bccode,bcvalue,sk,f);    
% 
%
  usparse = (sk\f);
%
clear sk;
clear f;


  u = full(usparse);
  
  clear usparse;
%
%remove tiny displacements caused by blasting  (cosmetic)
%
tinydispl = 1.0e-7*norm(nomdisplacement);
% for i=1:3*numnp
%     if(abs(u(i))<tinydispl) 
%         u(i) = 0;
%     end
% end

% Post-process to compute strain, stress, and mechanical metrics
%
%  first compute some geometric quantities for grains (just grain volume
%  for  now)
%
% grain_volumes = zeros(1,ngrains);
% grain_numels = zeros(1,ngrains);
% 
% for   iele =1:1:numel
%          
% igrain = grains(iele);
% 
% grain_volumes(igrain) = grain_volumes(igrain)+element_volumes(iele);
% grain_numels(igrain)  = grain_numels(igrain)+1;
% 
% end
% 
% sample_volume = sum(grain_volumes);


%   now compute the strain, stress and deviatoric stress
%
epsilon = zeros(6,nqptv,numel);
epsilon_ave= zeros(6,numel);
epsilon_average = zeros(6,1);

for   iele =1:1:numel
         
igrain = grains(iele);
%    
%compute strain and stress at the element quadrature points
%
[elepsilon_ave]=PostProcessForStrains(iele,nnpe,nqptv,wtqp,np,x,y,z,sfac,dndxi,dndet,dndze,u);

%
% compute running sum for the average stress 
%
epsilon_average(:) = epsilon_average(:) + elepsilon_ave(:)*element_volumes(iele);
%
%save elemental quantities in global arrays
%
epsilon_ave(:,iele) = elepsilon_ave(:);

end

 
 




