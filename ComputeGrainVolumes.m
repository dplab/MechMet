function [grain_volumes,grain_weights,grain_numels,sample_volume] = ComputeGrainVolumes(numel,np,x,y,z,ngrains,grains,phases,shafac_structure)
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

element_volumes = zeros(1,numel);
%
%
% Compute the elemental stiffness and force arrays and assemble
%
for   iele =1:1:numel
  
  [dndx,dndy,dndz,detj,elvolume] = SFDerivatives(iele,nnpe,nqptv,wtqp,np,x,y,z,dndxi,dndet,dndze);
  
  element_volumes(iele) = elvolume;

end

%
grain_volumes = zeros(1,ngrains);
grain_numels = zeros(1,ngrains);

for   iele =1:1:numel
         
igrain = grains(iele);

grain_volumes(igrain) = grain_volumes(igrain)+element_volumes(iele);
grain_numels(igrain)  = grain_numels(igrain)+1;

end

sample_volume = sum(grain_volumes);
grain_weights = grain_volumes/sample_volume;


 
 




