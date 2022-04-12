function [sx_compliance] =   Stiffness2Compliance(crystal_type,sx_moduli)

% compute the compliance components from the stiffness components 
%   input:  crystal type and moduli as designated below
%   output:  single crystal compliance components

% two symmetries are now supported: cubic and hexagonal
% In the s-x moduli file, specify  3 for cubic and 6 for hexagonal 
%  (moduli needed:
%    cubic:   c_11, c_12, and c_44
%    hexagonal:  c_11, c_12, c_13, and c_44 (c_33 and c_66 are computed))
%
%    Use strength of materials convention for shears:
%      sig_ij = c_44*gamma_ij  =  c_44* (2*eps_ij)  (i not equal j)
%



if(crystal_type == 3)
    
c11 = sx_moduli(1);
c12 = sx_moduli(2);
c44 = sx_moduli(3);

s11 = (c11+c12)/(c11*c11+c11*c12-2*c12*c12);
s12 = -c12/(c11*c11+c11*c12-2*c12*c12);
s44 = 1/c44;

sx_compliance = [s11 s12 s44];

elseif(crystal_type == 4)
    
c11 = sx_moduli(1);
c12 = sx_moduli(2);
c44 = sx_moduli(3);

s11 = (c11+c12)/(c11*c11+c11*c12-2*c12*c12);
s12 = -c12/(c11*c11+c11*c12-2*c12*c12);
s44 = 1/c44;

sx_compliance = [s11 s12 s44];
          
elseif(crystal_type==6)

c11 = sx_moduli(1);
c12 = sx_moduli(2);
c13 = sx_moduli(3);
c44 = sx_moduli(4);

c66 = (c11-c12)/2;
c33 =  c11+c12-c13;

r_crystal =  [c11 c12 c13  0    0    0; ...
                 c12 c11 c13  0    0    0; ...
                 c13 c13 c33  0    0    0; ...
                  0   0    0  c44  0    0; ...
                  0   0    0  0    c44  0; ...
                  0   0    0  0    0    c66];
              
s_crystal = inv(r_crystal);
             
s11 = s_crystal(1,1);
s12 = s_crystal(1,2);
s13 = s_crystal(1,3);
s44 = s_crystal(4,4);

sx_compliance = [s11 s12 s13 s44];            
end
%
