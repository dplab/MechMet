function [sx_moduli,residual] = MechModuli()
%  program MechModuli
%
%  Sister code to MechMet 
%
%  Purpose is to determine values of elastic moduli from HEXD lattice
%  strain data.  Using optimization method, this code varies values of elastic moduli
%  to find set of moduli that provides best match to measured lattice strains in selected grains
%
%  Uses same finite element methodology as MechMet to determine elasticity solution.  
%  Uses a univariant descent method with line search for the optimation.
%
%    Input:   (1a)  Neper-generated grain, orientation and mesh files for V3
%          or (1b)  Neper-generated msh for V4 (recommended choice)
%         and (2)  Lattice strain CHANGES (6xnumber of grains; SOM convention for shears) and corresponding macroscopic load CHANGE
%         and (3)  Binary flag to define the set of grains to be used in moduli determination (optional)
%                 (row matrix (1xnumber of grains) with values of zero or unity.  Used if the weight flag is set to zero.  
%         and (4)  the type of crystal symmetry (3 for cubic fcc, 4 for cubic bcc, or 6 for hexagonal (not yet verified)
%         and (5)  an estimate of the anisotropic ratio used in defining the initial guess for the set of moduli
%         and (6) the type of loading (x-, y- or z-direction tensile loading are the expected cases).
%
%    Output:  Set of single cyrstal elastic moduli and plots of the convergence history
%          

status = ' Starting procedure:  Read files with input data.'


neperversionnumber = input('Neper Version: Enter 3 for V3 or 4 for V4: ', 's');
neperversionnumber = str2num(neperversionnumber);

if(neperversionnumber==3) 
%  Read Neper phase, grain and orientation data

nepermicroname = input('Basename of microstructure files (with extensions .grain and .kocks): ', 's');

[grains, phases, orientations] = ReadNeperMicrostructure(nepermicroname);
bungeangles = BungeOfKocks(orientations,'degrees');
rotations = RMatOfBunge(bungeangles,'degrees');
ngrains=max(grains);
clear bungeangles;

microstructurefile = struct('phases',phases,'grains',grains,'orientations',orientations,'rotations',rotations);

% Read Neper mesh files
nepermeshname = input('Basename (with extension .mesh) of mesh file: ', 's');

meshfile = ReadNeperMesh(nepermeshname);

elseif(neperversionnumber==4)

nepermshname = input('Basename (with extension .msh) of the combined mesh  and microstructure file: ', 's');

[meshfile,microstructurefile] = ReadNeperMshFile(nepermshname);

phases = microstructurefile.phases;
grains = microstructurefile.grains;
rotations = microstructurefile.rotations;
orientations = microstructurefile.orientations;
ngrains=max(grains);

    
end

%
% Set mesh size variables from mesh coordinate and connectivity arrays
% 

[m n]   = size(meshfile.crd);
numeq = 3*n;
numnp = n;

[m n]   = size(meshfile.con);
numel = n;

np = meshfile.con';

%meshfile.crd = meshfile.crd*1e-6;


%
%define local variables for coordinates
%
x = meshfile.crd(1,:)';
y = meshfile.crd(2,:)';
z = meshfile.crd(3,:)';

%
%  set up reference array that gives phase of each grain
%

phaseofgrain=zeros(1,ngrains);

for iele=1:numel
 phasenum = phases(iele);
 igrain = grains(iele);
 phaseofgrain(igrain) = phasenum;
end

%
% Set up the finite element shape function arrays and store structures
%

shafac_structure = ShapeFunctionArrays();
shafac_surf_structure = ShapeFunctionArrays_Surface();


[grain_volumes,grain_weights,grain_numels,sample_volume] = ComputeGrainVolumes(numel,np,x,y,z,ngrains,grains,phases,shafac_structure);


%  
% weight_flag gives option for weighting grain errors in the residual:
%  weight_flag = 0 uses grainlist to mask (to use or ignore) grains 
%   (Use 0 for mask ; use 1 for relative grain volumes, which are computed within); 
%   Could be modified to include a relative weight based on confidence of data.
%

weightchoice = input('Input value of weight_flag: ', 's');

weight_flag = str2num(weightchoice);



if(weight_flag==1)
grain_wts = grain_weights;
else 

grainlistfilename = input('Basename of grainlist (.mat file with structure having .grainlist as a field): ', 's');  
grainliststructure = load([grainlistfilename, '.mat']);   
grainlist = grainliststructure.grainlist;    
numactivegrains = sum(grainlist);
grain_wts = grainlist/numactivegrains;
end


grainstrainsfilename = input('Basename of grain strains  (.mat file with structure having .grainstrains as a field): ', 's');  
grainstrainsstructure = load([grainstrainsfilename, '.mat']);   
grainstrains = grainstrainsstructure.grainstrains;  

grainstrainave_1 = grainstrains(1,:)*(grain_wts(:));
grainstrainave_2 = grainstrains(2,:)*(grain_wts(:));
grainstrainave_3 = grainstrains(3,:)*(grain_wts(:));


grainstrainave_4 = grainstrains(4,:)*(grain_wts(:));
grainstrainave_5 = grainstrains(5,:)*(grain_wts(:));
grainstrainave_6 = grainstrains(6,:)*(grain_wts(:));

% 
% 
% Set the loading type
%
message =  ['Loading code is defined as follows: '];
message1 =  ['loadcode = 1  --  x-tension'];
message2 =  ['loadcode = 2  --  y-tension'];
message3 =  ['loadcode = 3  --  z-tension'];
disp(message)
disp(message1)
disp(message2)
disp(message3)

loadcodes = input('Loadcode: ', 's');

loadcode = str2num(loadcodes);

loadvalue = input('Input the value of the load increment: ', 's');

load_increment = str2num(loadvalue);

        volumetric_strain = (grainstrainave_1+grainstrainave_2+grainstrainave_3);
switch loadcode
    case 1
        axial_strain = grainstrainave_1;
        lateral_strain = (grainstrainave_2+grainstrainave_3)/2;
    case 2
        axial_strain = grainstrainave_2;
        lateral_strain = (grainstrainave_1+grainstrainave_3)/2;
    case 3
        axial_strain = grainstrainave_3;
        lateral_strain = (grainstrainave_1+grainstrainave_2)/2;
end


% 
% Set up the boundary condition arrays for the specified loading mode
%


nominalstrain = 0.001;
lenxyz = [(max(x) - min(x)) (max(y) - min(y)) (max(z) - min(z))];
nomdisplacement = nominalstrain*lenxyz;

area_x = (max(y) - min(y))*(max(z) - min(z));
area_y = (max(x) - min(x))*(max(z) - min(z));
area_z = (max(x) - min(x))*(max(y) - min(y));
switch loadcode
    case 1 
        area = area_x;
    case 2 
        area = area_y;
    case 3 
        area = area_z;
    
end
traction_increment = load_increment/area;
%

[bccode,bcvalue,nomstress] = SetBCs_Moduli(meshfile,loadcode,nomdisplacement,neperversionnumber);


[MFsr,MFsv,Mfc] = ComputeForceMatrices(meshfile,shafac_surf_structure,numnp,loadcode, traction_increment);


%
% Input the crystal-type .  
%  
%   Do this for each of the phases in the neper file input

% two symmetries are now supported: cubic and hexagonal
% In the s-x moduli file, specify  3 (fcc) or 4 (bcc) for cubic and 6 for hexagonal 

nphases = max(phases);

message = ['Expecting material parameters for   ', num2str(nphases), ' phase(s)'];
disp(message)

for iphase = 1:nphases
    
message = ['Working on Phase ', num2str(iphase)];
disp(message)
    
ctype_string = input('Input crystal type: ', 's');

crystal_type = str2num(ctype_string);


A_string = input('Input estimate of anisotropic ratio, A: ', 's');

A_ratio = str2num(A_string);

%
%evaulate the single crystal stiffness and slip system strengthl
%                                                    
[sx_moduli,r_matrix_initial] =   SX_Moduli_InitialGuess(crystal_type,axial_strain,lateral_strain,volumetric_strain,traction_increment,A_ratio);

message = ['Initial moduli estimated from average strains: '];
disp(message)
sx_moduli

crystal_type_all(iphase) = crystal_type;
r_matrix_all(:,:,iphase) = r_matrix_initial(:,:);

end

%
%compute and store rotated stiffness matrices for each grain
%
grain_rmatrix = zeros(6,6,ngrains);

for igrain = 1:ngrains

myrot = rotations(:,:,igrain);
phasenum = phaseofgrain(igrain);
r_matrix(:,:) = r_matrix_all(:,:,phasenum);

grain_rmatrix(:,:,igrain) =  RotateStiffnessMatrix(r_matrix,myrot);

end


for iphase = 1:nphases

crystal_type = crystal_type_all(iphase);
                                                  
[r_matrix] =   SX_Moduli_Iterations(crystal_type,sx_moduli);

r_matrix_all(:,:,iphase) = r_matrix(:,:);

end

grain_rmatrix = zeros(6,6,ngrains);

for igrain = 1:ngrains

myrot = rotations(:,:,igrain);
phasenum = phaseofgrain(igrain);
r_matrix(:,:) = r_matrix_all(:,:,phasenum);

grain_rmatrix(:,:,igrain) =  RotateStiffnessMatrix(r_matrix,myrot);

end


Fsr = MFsr;
Fsv = MFsv;
fc = Mfc;

[epsilon_ave] = ComputeFELatticeStrains(numel,numnp,np,x,y,z,ngrains,grains,phases,shafac_structure,grain_rmatrix,bccode,bcvalue,Fsr,Fsv,fc,nomdisplacement);

computedstrain = zeros(6,ngrains);
strainerror = zeros(6,ngrains);
for iele = 1:1:numel
    grnnum = grains(iele);
    computedstrain(:,grnnum) = computedstrain(:,grnnum) + epsilon_ave(:,iele);
end
for igrn=1:1:ngrains
    computedstrain(:,igrn) = computedstrain(:,igrn)/grain_numels(igrn);
end

strainerror = grainstrains - computedstrain;
for icomp =1:1:6
 strainerror(icomp,:) = strainerror(icomp,:).*grain_wts(:)';   
end
residual_vector = reshape(strainerror,[6*ngrains,1]);
residual = residual_vector'*residual_vector;
residual_last = residual;


message = ['Residual following initial solution using initial moduli: ', num2str(residual)];
disp(message)


computedstrainave_1 = computedstrain(1,:)*grain_wts(:);
computedstrainave_2 = computedstrain(2,:)*grain_wts(:);
computedstrainave_3 = computedstrain(3,:)*grain_wts(:);
% 


switch loadcode
    case 1
        axial_strain_comp = computedstrainave_1;
        lateral_strain_comp = (computedstrainave_2+computedstrainave_3)/2;
    case 2
        axial_strain_comp = computedstrainave_2;
        lateral_strain_comp = (computedstrainave_1+computedstrainave_3)/2;
    case 3
        axial_strain_comp = computedstrainave_3;
        lateral_strain_comp = (computedstrainave_1+computedstrainave_2)/2;
end

scalefactor = axial_strain/axial_strain_comp;

sx_moduli = sx_moduli/scalefactor;  

message = ['Rescaled moduli estimated from ratio of computed to average strains: '];
disp(message)
sx_moduli
% 
% now recompute the residual for the new set of moduli
%
                for iphase = 1:nphases

                crystal_type = crystal_type_all(iphase);

                [r_matrix_scaled] =   SX_Moduli_Iterations(crystal_type,sx_moduli);

                r_matrix_all(:,:,iphase) = r_matrix_scaled(:,:);

                end

                grain_rmatrix = zeros(6,6,ngrains);

                for igrain = 1:ngrains

                myrot = rotations(:,:,igrain);
                phasenum = phaseofgrain(igrain);
                r_matrix(:,:) = r_matrix_all(:,:,phasenum);

                grain_rmatrix(:,:,igrain) =  RotateStiffnessMatrix(r_matrix,myrot);

                end


                Fsr = MFsr;
                Fsv = MFsv;
                fc = Mfc;

                [epsilon_ave] = ComputeFELatticeStrains(numel,numnp,np,x,y,z,ngrains,grains,phases,shafac_structure,grain_rmatrix,bccode,bcvalue,Fsr,Fsv,fc,nomdisplacement);

                computedstrain = zeros(6,ngrains);
                strainerror = zeros(6,ngrains);
                for iele = 1:1:numel
                    grnnum = grains(iele);
                    computedstrain(:,grnnum) = computedstrain(:,grnnum) + epsilon_ave(:,iele);
                end
                for igrn=1:1:ngrains
                    computedstrain(:,igrn) = computedstrain(:,igrn)/grain_numels(igrn);
                end

                strainerror = grainstrains - computedstrain;

                for icomp =1:1:6
                 strainerror(icomp,:) = strainerror(icomp,:).*grain_wts(:)';   
                end

                residual_vector = reshape(strainerror,[6*ngrains,1]);
                residual = residual_vector'*residual_vector;residual_last = residual;


message = ['Residual following rescaling using computed average strains: ', num2str(residual)];
disp(message)



%
% perform optimization to refine values of the moduli by reducting strain residual

  sizesx = size(sx_moduli);
  num_delta_sij = sizesx(2);
  
residual_target = 1.0e-5*residual
iteration_count = 0;
iteration_limit = 40;

maxchange = 0.2;
delta = 0.04;  

residual_history = zeros(iteration_limit*num_delta_sij,1);
jacobian_history = zeros(iteration_limit*num_delta_sij,1);
sx_moduli_history = zeros(iteration_limit*num_delta_sij,num_delta_sij);
sx_compliance_history = zeros(iteration_limit*num_delta_sij,num_delta_sij);
% 
 while(residual > residual_target)
%     
%     
     if(iteration_count> iteration_limit*num_delta_sij)
        save('residual_history.mat','residual_history');
        save('jacobian_history.mat','jacobian_history');
        save('sx_moduli_history.mat','sx_moduli_history');        
        save('sx_compliance_history.mat','sx_compliance_history');
        
        PlotModuliConvergence(residual_history,sx_moduli_history,sx_compliance_history)

        message =  ['Maximum iterations reached.  Hit Control C to exit. '];
        disp(message)
        pause
     end
% 

  for idelta = 1:1:num_delta_sij
      
      
  iteration_count = iteration_count+1;
  
    
  delta_ij = ones(1, num_delta_sij);
  delta_ij(idelta) = 1.0+delta;
      
  [sx_compliance] =   Stiffness2Compliance(crystal_type,sx_moduli);
 
  sx_compliance_delta = sx_compliance.*delta_ij;
  
  [sx_moduli_delta] =   Compliance2Stiffness(crystal_type,sx_compliance_delta);
  
  
for iphase = 1:nphases

crystal_type = crystal_type_all(iphase);
                                                  
[r_matrix_delta] =   SX_Moduli_Iterations(crystal_type,sx_moduli_delta);

r_matrix_all(:,:,iphase) = r_matrix_delta(:,:);

end

grain_rmatrix = zeros(6,6,ngrains);

for igrain = 1:ngrains

myrot = rotations(:,:,igrain);
phasenum = phaseofgrain(igrain);
r_matrix(:,:) = r_matrix_all(:,:,phasenum);

grain_rmatrix(:,:,igrain) =  RotateStiffnessMatrix(r_matrix,myrot);

end

 Fsr = MFsr;
 Fsv = MFsv;
 fc = Mfc;
% 
 [epsilon_ave] = ComputeFELatticeStrains(numel,numnp,np,x,y,z,ngrains,grains,phases,shafac_structure,grain_rmatrix,bccode,bcvalue,Fsr,Fsv,fc,nomdisplacement);
% 
% %
% 
delta_strain = zeros(6,ngrains);
for iele = 1:1:numel
    grnnum = grains(iele);
    delta_strain(:,grnnum) = delta_strain(:,grnnum) + epsilon_ave(:,iele);
end
for igrn=1:1:ngrains
    delta_strain(:,igrn) = delta_strain(:,igrn)/grain_numels(igrn);
end

deps_dsij = -(delta_strain - computedstrain)/(sx_compliance(idelta)*delta);

gradient_vector = reshape(deps_dsij,[6*ngrains,1]);

jacobian_sij = 2*residual_vector'*gradient_vector;

delta_sij_nr = -residual/jacobian_sij;

  
ilinecheck = 1;
  for ilinesearch = 1:1:12
      
    if(ilinecheck==1) 
      
            compliance_change = maxchange/(2^(ilinesearch-1)); 
            
            delta_sij_ls = compliance_change*sign(delta_sij_nr)*abs(sx_compliance(idelta));

            sx_compliance_delta(idelta) = sx_compliance(idelta) + delta_sij_ls;

            [sx_moduli] =   Compliance2Stiffness(crystal_type,sx_compliance_delta);


                % 
                % now recompute the residual for the new set of moduli
                for iphase = 1:nphases

                crystal_type = crystal_type_all(iphase);

                [r_matrix_scaled] =   SX_Moduli_Iterations(crystal_type,sx_moduli);

                r_matrix_all(:,:,iphase) = r_matrix_scaled(:,:);

                end

                grain_rmatrix = zeros(6,6,ngrains);

                for igrain = 1:ngrains

                myrot = rotations(:,:,igrain);
                phasenum = phaseofgrain(igrain);
                r_matrix(:,:) = r_matrix_all(:,:,phasenum);

                grain_rmatrix(:,:,igrain) =  RotateStiffnessMatrix(r_matrix,myrot);

                end


                Fsr = MFsr;
                Fsv = MFsv;
                fc = Mfc;

                [epsilon_ave] = ComputeFELatticeStrains(numel,numnp,np,x,y,z,ngrains,grains,phases,shafac_structure,grain_rmatrix,bccode,bcvalue,Fsr,Fsv,fc,nomdisplacement);

                computedstrain = zeros(6,ngrains);
                strainerror = zeros(6,ngrains);
                for iele = 1:1:numel
                    grnnum = grains(iele);
                    computedstrain(:,grnnum) = computedstrain(:,grnnum) + epsilon_ave(:,iele);
                end
                for igrn=1:1:ngrains
                    computedstrain(:,igrn) = computedstrain(:,igrn)/grain_numels(igrn);
                end

                strainerror = grainstrains - computedstrain;

                for icomp =1:1:6
                 strainerror(icomp,:) = strainerror(icomp,:).*grain_wts(:)';   
                end

                residual_vector = reshape(strainerror,[6*ngrains,1]);
                residual = residual_vector'*residual_vector;

     if(residual < residual_last)
       ilinecheck = 0;  
     end

      end
  end

residual_last = residual;
residual_history(iteration_count) = residual;
jacobian_history(iteration_count) = jacobian_sij;
sx_compliance = sx_compliance_delta;

sx_moduli_history(iteration_count,:) = sx_moduli(:);
sx_compliance_history(iteration_count,:) = sx_compliance(:);

  end
  
  
[sx_moduli] =   Compliance2Stiffness(crystal_type,sx_compliance);


message = ['Residual for iteration ', num2str(iteration_count),  ' is: ', num2str(residual)];
disp(message)
sx_moduli
sx_compliance

%

 end


        save('residual_history.mat','residual_history');
        save('jacobian_history.mat','jacobian_history');
        save('sx_moduli_history.mat','sx_moduli_history');
        save('sx_compliance_history.mat','sx_compliance_history');
        
        PlotModuliConvergence(residual_history,sx_moduli_history,sx_compliance_history)

end
