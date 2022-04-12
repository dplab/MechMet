function [sx_moduli] =   Compliance2Stiffness(crystal_type,sx_compliance)

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
    
s11 = sx_compliance(1);
s12 = sx_compliance(2);
s44 = sx_compliance(3);

c11 = (s11+s12)/((s11-s12)*(s11+2*s12));
c12 = -s12/((s11-s12)*(s11+2*s12));
c44 = 1/s44;

sx_moduli = [c11 c12 c44];

elseif(crystal_type == 4)
    
s11 = sx_compliance(1);
s12 = sx_compliance(2);
s44 = sx_compliance(3);

c11 = (s11+s12)/((s11-s12)*(s11+2*s12));
c12 = -s12/((s11-s12)*(s11+2*s12));
c44 = 1/s44;

sx_moduli = [c11 c12 c44];
          
elseif(crystal_type==6)

s11 = sx_compliance(1);
s12 = sx_compliance(2);
s13 = sx_compliance(3);
s44 = sx_compliance(4);

s66 = (s11-s12)/2;
s33 =  s11+s12-s13;

s_crystal =  [s11 s12 s13  0    0    0; ...
                 s12 s11 s13  0    0    0; ...
                 s13 s13 s33  0    0    0; ...
                  0   0    0  s44  0    0; ...
                  0   0    0  0    s44  0; ...
                  0   0    0  0    0    s66];
              
c_crystal = inv(s_crystal);
             
  
sx_compliance = [c11 c12 c13 c44];            
end
%
