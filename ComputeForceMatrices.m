
function [Fsr,Fsv,fc] = ComputeForceMatrices(meshfile,shafac_surf_structure,numnp,loadcode, traction_increment)
%
% computes the grain-average lattice strain over finite element mesh 
%          

ssfac = shafac_surf_structure.ssfac;
dnda = shafac_surf_structure.dnda;
dndb = shafac_surf_structure.dndb;
nqpts = shafac_surf_structure.nqpts;
swts = shafac_surf_structure.swts;
nnps = shafac_surf_structure.nnps;

Fsr=zeros(3*numnp,1); 
Fsv=zeros(3*numnp,1); 
fc=1; 

sft = zeros(3*numnp,1);
xlocal = zeros(6,1);
ylocal = zeros(6,1);
zlocal = zeros(6,1);

for isurf = 1:1:6
    
    nsurele = meshfile.nps_sizes(isurf);
    surf_trac = [0 0 0]';
    
    switch isurf
        case 1
            nps_n = meshfile.nps_1;
        case 2
            nps_n = meshfile.nps_2;
            if(loadcode==1) 
                surf_trac = [traction_increment 0 0 ]';
            end
        case 3
            ps_n = meshfile.nps_3;
        case 4
            nps_n = meshfile.nps_4;
            if(loadcode==2) 
                surf_trac = [ 0 traction_increment 0 ]';
            end
        case 5
            nps_n = meshfile.nps_6;
        case 6
            nps_n = meshfile.nps_6;
            if(loadcode==3) 
                surf_trac = [0 0 traction_increment]';
            end
     end
     
        
    if(norm(surf_trac) ~= 0 )     
        
        
      for jele = 1:1:nsurele
          
          solid_element_num = nps_n(jele,1);
          nps_l([1:6]) = nps_n(jele,[2:7]);
          nps([ 1 2 3 4 5 6]) = nps_l([6 1 4 2 5 3]);
          
          for ksnp = 1:1:6
          surfnpnum = nps(ksnp);    
          xlocal(ksnp) = meshfile.crd(1,surfnpnum);
          ylocal(ksnp) = meshfile.crd(2,surfnpnum);
          zlocal(ksnp) = meshfile.crd(3,surfnpnum);
          end   
    
        [n rjs] = SurfJacobian(jele,nnps,nqpts,dnda,dndb,nps,xlocal,ylocal,zlocal);
    
        bigNsurf = BigNSmat(nnps,nqpts,ssfac);
    
        fes = SurfForce(nnps,nqpts,swts,bigNsurf,rjs,surf_trac);
        
        [Fsr,Fsv,fc] = AssembleSurfTractions(jele,nnps,nps,fes,Fsr,Fsv,fc);

      end
      
    end 
    
end

end
