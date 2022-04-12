function [sx_moduli_initial,r_crystal] =   SX_Moduli_InitialGuess(crystal_type,axial_strain,lateral_strain,volumetric_strain,traction_increment,A)

% read in  initial guesses for single crystal constants and return the material stiffness 
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
options = optimoptions(@fmincon,'MaxIterations',10,'MaxFunctionEvaluations',1000000,'Display','iter')
options.MaxFunctionEvaluations = 1000000;
options.OptimalityTolerance = 1.0e-12;
options.StepTolerance  = 1.0e-6;

% v=v*norm_data;


if(crystal_type == 3)

      YoungMod = traction_increment/axial_strain;
      PoissonRatio = -lateral_strain/axial_strain;
      
      ci(1) = YoungMod*(1-PoissonRatio)/(1-PoissonRatio-2*PoissonRatio^2);
      ci(2) = YoungMod*PoissonRatio/(1-PoissonRatio-2*PoissonRatio^2);
      ci(3) = (ci(1)-ci(2))*A/2;
      
       E = traction_increment/axial_strain;
       B= traction_increment/(3.0*volumetric_strain);
       
       Aeq = [ 1 2 0; A -A -2 ];
       Beq = [B ; 0 ];

      c =zeros(3,1);              

     [c] = fmincon(@(c)((c(1)-c(2)+3*c(3))*(c(1)+2*c(2))-E*(2*c(1)+3*c(2)+c(3)))^2,ci,[],[],Aeq,Beq,[],[],[],options);
       
      c11 = c(1);
      c12 = c(2);
      c44 = c(3);

      
      
      sx_moduli_initial = [c11 c12 c44];

r_crystal =  [c11 c12 c12  0    0    0; ...
              c12 c11 c12  0    0    0; ...
              c12 c12 c11  0    0    0; ...
              0   0    0  c44  0    0; ...
              0   0    0  0    c44  0; ...
              0   0    0  0    0    c44];
          

elseif(crystal_type == 4)

      YoungMod = traction_increment/axial_strain;
      PoissonRatio = -lateral_strain/axial_strain;
      
      ci(1) = YoungMod*(1-PoissonRatio)/(1-PoissonRatio-2*PoissonRatio^2);
      ci(2) = YoungMod*PoissonRatio/(1-PoissonRatio-2*PoissonRatio^2);
      ci(3) = (ci(1)-ci(2))*A/2;
      
       E = traction_increment/axial_strain;
       B= traction_increment/(3.0*volumetric_strain);
       
       Aeq = [ 1 2 0; A -A -2 ];
       Beq = [B ; 0 ];

      c =zeros(3,1);              

     [c] = fmincon(@(c)((c(1)-c(2)+3*c(3))*(c(1)+2*c(2))-E*(2*c(1)+3*c(2)+c(3)))^2,ci,[],[],Aeq,Beq,[],[],[],options);
       
      c11 = c(1);
      c12 = c(2);
      c44 = c(3);
      sx_moduli_initial = [c11 c12 c44];

r_crystal =  [c11 c12 c12  0    0    0; ...
              c12 c11 c12  0    0    0; ...
              c12 c12 c11  0    0    0; ...
              0   0    0  c44  0    0; ...
              0   0    0  0    c44  0; ...
              0   0    0  0    0    c44];
          
          

          
elseif(crystal_type==6)

      A=1;
      YoungMod = traction_increment/axial_strain;
      PoissonRatio = -lateral_strain/axial_strain;
      c11 = YoungMod*(1-PoissonRatio)/(1-PoissonRatio-2*PoissonRatio^2);
      c12 = YoungMod*PoissonRatio/(1-PoissonRatio-2*PoissonRatio^2);
      c13=c12;
      c44 = (c11-c12)*A/2;
      
      sx_moduli_initial = [c11 c12 c13 c44];
      
c66 = (c11-c12)/2;
c33 =  c11+c12-c13;

r_crystal =  [c11 c12 c13  0    0    0; ...
                 c12 c11 c13  0    0    0; ...
                 c13 c13 c33  0    0    0; ...
                  0   0    0  c44  0    0; ...
                  0   0    0  0    c44  0; ...
                  0   0    0  0    0    c66];
             
              
end
%
end
function f = voigtestimate(c,B,E,A)
 f3 = [c(1)+2*c(2)-B; 2*c(3)-A*(c(1)-c(2)); (c(1)-c(2)+3*c(3))*(c(1)+2*c(2))-E*(2*c(1)+3*c(2)+c(3))];
 f = f3'*f3;
end
     function f = voigtestimate2(c,E)
      f = ((c(1)-c(2)+3*c(3))*(c(1)+2*c(2))-E*(2*c(1)+3*c(2)+c(3)))^2;
     end
