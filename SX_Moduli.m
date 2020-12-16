function [crystal_type,r_crystal,sss,sssr,covera] =   SX_Moduli(sx_filename)
% read in single crystal constants and return the material stiffness 
%  matrix in crystal reference frame
%   input:  crystal type and moduli as designated below
%   output:  single crystal material stiffness matrix

% two symmetries are now supported: cubic and hexagonal
% In the s-x moduli file, specify  3 for cubic and 6 for hexagonal 
%  (moduli needed:
%    cubic:   c_11, c_12, and c_44
%    hexagonal:  c_11, c_12, c_13, and c_44 (c_33 and c_66 are computed))
%
%    Use strength of materials convention for shears:
%      sig_ij = c_44*gamma_ij  =  c_44* (2*eps_ij)  (i not equal j)
%
% Slip system strengths:
%    For FCC, enter the strength of the slip systems (one value
%    for all) and slip system type (1 for fcc; 2 for bcc).
%    For BCC, enter 4 values: Strength of {110}, {112}, {123} 
%    and {134} slip planes.
%    For HCP, enter 4 values:  Strength of basal, strength of prismatic,
%    strength of pyramidal, and c over a ratio of the unit cell.

modulifile = [sx_filename,'.matl'];
fid = fopen(modulifile,'r');
    
crystal_type = fscanf(fid, '%d', 1);

if(crystal_type == 3)
parms = fscanf(fid, '%f %f %f', 3);
c11 = parms(1);
c12 = parms(2);
c44 = parms(3);

r_crystal =  [c11 c12 c12  0    0    0; ...
              c12 c11 c12  0    0    0; ...
              c12 c12 c11  0    0    0; ...
              0   0    0  c44  0    0; ...
              0   0    0  0    c44  0; ...
              0   0    0  0    0    c44];
          
          

ssval = fscanf(fid, '%f %f', 2);
ssst = [
    ssval(1) ssval(1) ssval(1) ssval(1) ssval(1) ssval(1) ... 
    ssval(1) ssval(1) ssval(1) ssval(1) ssval(1) ssval(1) ... 
    ssval(1) ssval(1) ssval(1) ssval(1) ssval(1) ssval(1) 
    ];
sss = ssst';
sssr = sss/ssval(1);
covera = ssval(2);

elseif(crystal_type == 4)
parms = fscanf(fid, '%f %f %f', 3);
c11 = parms(1);
c12 = parms(2);
c44 = parms(3);

r_crystal =  [c11 c12 c12  0    0    0; ...
              c12 c11 c12  0    0    0; ...
              c12 c12 c11  0    0    0; ...
              0   0    0  c44  0    0; ...
              0   0    0  0    c44  0; ...
              0   0    0  0    0    c44];
          
          

ssval = fscanf(fid, '%f %f %f %f %f', 5);
ssst = [
    ssval(1) ssval(1) ssval(1) ssval(1) ssval(1) ssval(1) ... 
    ssval(1) ssval(1) ssval(1) ssval(1) ssval(1) ssval(1) ... 
    ssval(2) ssval(2) ssval(2) ssval(2) ssval(2) ssval(2) ...
    ssval(2) ssval(2) ssval(2) ssval(2) ssval(2) ssval(2) ...
    ssval(3) ssval(3) ssval(3) ssval(3) ssval(3) ssval(3) ...
    ssval(3) ssval(3) ssval(3) ssval(3) ssval(3) ssval(3) ...
    ssval(3) ssval(3) ssval(3) ssval(3) ssval(3) ssval(3) ...
    ssval(3) ssval(3) ssval(3) ssval(3) ssval(3) ssval(3) ...
    ssval(4) ssval(4) ssval(4) ssval(4) ssval(4) ssval(4) ...
    ssval(4) ssval(4) ssval(4) ssval(4) ssval(4) ssval(4) ...
    ssval(4) ssval(4) ssval(4) ssval(4) ssval(4) ssval(4) ...
    ssval(4) ssval(4) ssval(4) ssval(4) ssval(4) ssval(4) ...
    ];
sss = ssst';
sssr = sss/ssval(1);
covera = ssval(5);
          
elseif(crystal_type==6)

parms = fscanf(fid, '%f %f %f %f', 4);
c11 = parms(1);
c12 = parms(2);
c13 = parms(3);
c44 = parms(4);
c66 = (c11-c12)/2;
c33 =  c11+c12-c13;

r_crystal =  [c11 c12 c13  0    0    0; ...
                 c12 c11 c13  0    0    0; ...
                 c13 c13 c33  0    0    0; ...
                  0   0    0  c44  0    0; ...
                  0   0    0  0    c44  0; ...
                  0   0    0  0    0    c66];
              
                        
ssval = fscanf(fid, '%f %f %f %f', 4);
ssst = [ssval(1) ssval(1) ssval(1) ssval(2) ssval(2) ssval(2) ...
    ssval(3) ssval(3) ssval(3) ssval(3) ssval(3) ssval(3) ...
    ssval(3) ssval(3) ssval(3) ssval(3) ssval(3) ssval(3)
    ];
covera = ssval(4);
sss = ssst';
sssr = sss/ssval(1);
              
end
%
status = fclose(fid);
