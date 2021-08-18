function [ReturnStatus] = MechMet()
%  program MechMet
%
%  Compiled by Paul Dawson from codes and scripts from MAE courses, OdfPF, 
%   pre- and post-processing tools for FEpX, other research projects
%   in the Deformation Processes Lab at Cornell University.
%
%  Driver program for 3-d elasticity solution of virtual polycrystals 
%    Intent:   provide linear elastic solution for polycrystal samples using mesh
%       files generated by Neper (especially from HEDM data)
%
%    Input:   (1)  Neper-generated grain, orientation and mesh files for V3
%          or (2)  Neper-generated msh for V4
%         and (3)  User-generated crystal property files
%
%    Output:  Using anisotropic elasticity field equations:
%               (1) displacement field for designated loading mode
%               (2) elastic strain field
%               (3) stress field 
%               (4) grain-smoothed strain and stress fields
%               (5) mechanical metrics (eg, RE2SX ratio, Y2E, SchmidFactor)
%          

ReturnStatus = 'Failed';

status = ' Starting simulation:  Read files with input data.'


neperversionnumber = input('Neper Version: Enter 3 for V3 or 4 for V4: ', 's');
neperversionnumber = str2num(neperversionnumber);

if(neperversionnumber==3) 
%  Read Neper phase, grain and orientation data

nepermicroname = input('Base name of microstructure files (with extensions .grain and .kocks): ', 's');

[grains, phases, orientations] = ReadNeperMicrostructure(nepermicroname);
bungeangles = BungeOfKocks(orientations,'degrees');
rotations = RMatOfBunge(bungeangles,'degrees');
ngrains=max(grains);
clear bungeangles;

microstructurefile = struct('phases',phases,'grains',grains,'orientations',orientations,'rotations',rotations);

% Read Neper mesh files
nepermeshname = input('Base name (with extension .mesh) of mesh file: ', 's');

meshfile = ReadNeperMesh(nepermeshname);

elseif(neperversionnumber==4)

nepermshname = input('Base name (with extension .msh) of the combined mesh  and microstructure file: ', 's');

[meshfile,microstructurefile] = ReadNeperMshFile(nepermshname);

phases = microstructurefile.phases;
grains = microstructurefile.grains;
rotations = microstructurefile.rotations;
orientations = microstructurefile.orientations;
ngrains=max(grains);

    
end


% Set mesh size variables from mesh coordinate and connectivity arrays
[m n]   = size(meshfile.crd);
numeq = 3*n;
numnp = n;

[m n]   = size(meshfile.con);
numel = n;

np = meshfile.con';
coords = meshfile.crd;



%
%define local variables for coordinates
%

x = coords(1,:)';
y = coords(2,:)';
z = coords(3,:)';

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
% Set up the finite element shape function arrays and store in structure
%

shafac_structure = ShapeFunctionArrays();

sfac = shafac_structure.sfac;
dndxi = shafac_structure.dndxi;
dndet = shafac_structure.dndet;
dndze = shafac_structure.dndze;
nqptv = shafac_structure.nqpt;
wtqp = shafac_structure.wtqp;
nnpe = shafac_structure.nnpe;
%


% Input the crystal-type, elastic moduli and slip systems strengths.  
%  
%   Do this for each of the phases in the neper file input

% two symmetries are now supported: cubic and hexagonal
% In the s-x moduli file, specify  3 for cubic and 6 for hexagonal 
%  (moduli needed
%    cubic:   c_11, c_12, and c_44
%    hexagonal:  c_11, c_12, c_13, and c_44 (c_33 and c_66 are computed))
%
% Slip system strengths:
%    For FCC, enter the strength of the slip systems (one value
%    for all) and slip system type (1 for fcc; 2 for bcc).
%    For BCC, enter 4 values: Strength of {110}, {112}, {123} 
%    and {134} slip planes.
%    For HCP, enter 4 values:  Strength of basal, strength of prismatic,
%    strength of pyramidal, and c over a ratio of the unit cell.
%
%   maxss should be increased if more that 18 slip systems are used for any
%   phase  (currently the number for hexagonal crystal type)
%   12 ss are used for cubics; slip systems 13-18 are zeros. 

%   Update: maxss was increased due to the 72 unique slip systems for BCC metals
%   when considering different slip strength on
%   12 {110}, 12 {112}, 24 {123}, 24 {134} systems.

nphases = max(phases);
maxss = 72;
crytal_type_all = zeros(1,nphases);
r_matrix_all = zeros(6,6,nphases);
sss_all = zeros(maxss,nphases);
sssr_all = zeros(maxss,nphases);
covera_all = ones(1,nphases);

message = ['Expecting material parameters for   ', num2str(nphases), ' phase(s)'];
disp(message)

for iphase = 1:nphases
    
message = ['Working on Phase ', num2str(iphase)];
disp(message)
    
sx_filename = input('Base name of single-crystal material file (with extension .matl): ', 's');
%
%evaulate the single crystal stiffness and slip system strengthl
%                                                    
[crystal_type,r_matrix,sss,sssr,covera] =   SX_Moduli(sx_filename,maxss);

crystal_type_all(iphase) = crystal_type;
r_matrix_all(:,:,iphase) = r_matrix(:,:);
sss_all(:,iphase) = sss(:);
sssr_all(:,iphase) = sssr(:);
covera_all(iphase) = covera;

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


% 
% Set the loading type
%
message =  ['Loading code is defined as follows: '];
message1 =  ['loadcode = 1 --  x-tension'];
message2 =  ['loadcode = 2  --  y-tension'];
message3 =  ['loadcode = 3  --  z-tension'];
message4 =  ['loadcode = 4  --  x-y biaxial tension'];
message5 =  ['loadcode = 5  --  y-z biaxial tension'];
message6 =  ['loadcode = 6  --  z-x biaxial tension'];
disp(message)
disp(message1)
disp(message2)
disp(message3)
disp(message4)
disp(message5)
disp(message6)

loadcodes = input('Loadcode: ', 's');

% the code is as follows:
%   loadcode = 1  --  x-tension
%   loadcode = 2  --  y-tension
%   loadcode = 3  --  z-tension
%   loadcode = 4  --  x-y biaxial tension
%   loadcode = 5  --  y-z biaxial tension
%   loadcode = 6  --  z-x biaxial tension

loadcode = str2num(loadcodes);


%
% allocate arrays for the global stiffness and force matrices
%  using sparse matrix solver 


Fsr=zeros(numel*3*nnpe,1); 
Fsv=zeros(numel*3*nnpe,1); 
Ksr=zeros(numel*9*nnpe*nnpe,1); 
Ksc=zeros(numel*9*nnpe*nnpe,1); 
Ksv=zeros(numel*9*nnpe*nnpe,1);
fc=1; 
kc=1;

%
%defines body forces to be zero 
%(potentially could be used to include residual stresses)
%
gamvec = [ 0 0 0 ];
element_volumes = zeros(1,numel);
%

status = ' Starting Element Integration and Assembly'
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


status = ' Finished Assembly '

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


status = ' Starting Displacement Solution'


u = zeros(numeq,1);

% 
% Set up the boundary condition arrays for the specified loading mode
%
%
%compute sample dimensions and set compute displacements to produce a
%nominal strain of 0.001
%
nominalstrain = 0.001;
lenxyz = [(max(x) - min(x)) (max(y) - min(y)) (max(z) - min(z))];
nomdisplacement = nominalstrain*lenxyz;
%
%set nominal strain for strength-to-stiffness computation
Eeff= nominalstrain;
%

[bccode,bcvalue,nomstress] = SetBCs(meshfile,loadcode,nomdisplacement,neperversionnumber);

%
% Apply boundary conditions
% 
  [sk,f] = BCSparse(bccode,bcvalue,sk,f);    
% 
%
  u=full(sk\f);
%
clear sk;
clear f;
%
%remove tiny displacements caused by blasting  (cosmetic)
%
tinydispl = 1.0e-7*norm(nomdisplacement);
for i=1:numeq
    if(abs(u(i))<tinydispl) 
        u(i) = 0;
    end
end

status = ' Finished Displacement Solution  --  Starting Post-Processing'

%
% Post-process to compute strain, stress, and mechanical metrics
%
%  first compute some geometric quantities for grains (just grain volume
%  for  now)
%
grain_volumes = zeros(1,ngrains);
grain_numels = zeros(1,ngrains);

for   iele =1:1:numel
         
igrain = grains(iele);

grain_volumes(igrain) = grain_volumes(igrain)+element_volumes(iele);
grain_numels(igrain)  = grain_numels(igrain)+1;

end

sample_volume = sum(grain_volumes);


%   now compute the strain, stress and deviatoric stress
%
epsilon = zeros(6,nqptv,numel);
sigma = zeros(6,nqptv,numel);
epsilon_ave= zeros(6,numel);
sigma_ave=zeros(6,numel);
sigma_average = zeros(6,1);
epsilon_average = zeros(6,1);

sigmaxs = zeros(6,1);
sigmins = zeros(6,1);
sigrangeqp = zeros(numel,1);

for   iele =1:1:numel
         
igrain = grains(iele);
rotated_rmatrix = grain_rmatrix(:,:,igrain);
%    
%compute strain and stress at the element quadrature points
%
[elepsilon,elepsilon_ave,elsigma,elsigma_ave,eldsigma]=PostProcess(iele,rotated_rmatrix,nnpe,nqptv,wtqp,np,x,y,z,sfac,dndxi,dndet,dndze,u);
%
% compute running sum for the average stress 
%
sigma_average(:) = sigma_average(:) + elsigma_ave(:)*element_volumes(iele);
epsilon_average(:) = epsilon_average(:) + elepsilon_ave(:)*element_volumes(iele);
%
%save elemental quantities in global arrays
%   (strain saved in elasticity notation)
%
for iqpt = 1:1:nqptv
sigma(:,iqpt,iele) = elsigma(:,iqpt);
epsilon_sm(1:6) = [elepsilon(1,iqpt);elepsilon(2,iqpt);elepsilon(3,iqpt);elepsilon(4,iqpt)/2.;elepsilon(5,iqpt)/2.;elepsilon(6,iqpt)/2.];
epsilon(:,iqpt,iele) = epsilon_sm(:);
end

for icomp =1:1:6
sigmaxs(icomp) = max(elsigma(icomp,:));
sigmins(icomp) = min(elsigma(icomp,:));
end
sigrangeqp(iele) = max(sigmaxs(:)-sigmins(:));

sigma_ave(:,iele) = elsigma_ave(:);
epsilon_ave(:,iele) = elepsilon_ave(:);

end

sigma_average = sigma_average/sample_volume;

[maxvalue,maxcomp] = max(sigma_average);


%
% begin smoothing of derivative quantities:  stress and strain 
% 
%       This uses a consistent penalty method applied grain-by-grain%

%   First, reconstuct the mesh without intragrain continuity (nodes on the
%       grain boundaries have separate numbers for each grain.
%



status = ' Re-numbering mesh for grain-by-grain stress and strain continuity. '



[newnumnp,newnp,newnpinv,newcoords,old2new,grain4np,meshfilegbg]=GrainByGrainMesh(numel,numnp,np,nnpe,coords,grains);



%
%
% Compute the elemental Consistent Penalty RHS and LHS arrays
%
% Average at nodal points
%



status = ' Starting Element-by-Element collocation of stress and strain, followed by averaging at nodes'



damping_factors = zeros(numel,1);
dampflagel = zeros(numel,1);
sigrangenp = zeros(numel,1);

sigmanpave = zeros(6,newnumnp);
epsilonnpave = zeros(6,newnumnp);
numelepernp = zeros(newnumnp,1);


numfields = 12;
feall = zeros(nnpe,numfields);

xnew = newcoords(1,:)';
ynew = newcoords(2,:)';
znew = newcoords(3,:)';

for   iele =1:1:numel
    
    
  sigmanp = zeros(6,nnpe);
  epsilonnp = zeros(6,nnpe);
  

    
dampflag = 0;
damp_iter =0;



while (dampflag ==0)
    
 damp_iter = damp_iter+1;
      
  igrain = grains(iele);
  dampfac = damping_factors(iele);
  
  [dndx,dndy,dndz,detj,elvolume] = SFDerivatives(iele,nnpe,nqptv,wtqp,newnp,xnew,ynew,znew,dndxi,dndet,dndze);
  
  bigN = sfac;

  bigNprime = BigNprimeMat(nnpe,nqptv,dndx,dndy,dndz);
  
  
  me = ElementLHS4ConPen(nnpe,nqptv,wtqp,bigN,bigNprime,detj,dampfac);

    
    for ifield = 1:1:numfields
    
       switch ifield
       case 1
            data = sigma(1,:,iele);
       case 2         
            data = sigma(2,:,iele);
       case 3  
            data = sigma(3,:,iele);
       case 4
            data = sigma(4,:,iele);  
       case 5
            data = sigma(5,:,iele);   
       case 6
            data = sigma(6,:,iele);
       case 7
            data = epsilon(1,:,iele);
       case 8         
            data = epsilon(2,:,iele);
       case 9  
            data = epsilon(3,:,iele);
       case 10
            data = epsilon(4,:,iele);  
       case 11
            data = epsilon(5,:,iele);   
       case 12
            data = epsilon(6,:,iele);
       end
       
    fe = ElementRHS4ConPen(nnpe,nqptv,wtqp,bigN,detj,data);
    
    feall(:,ifield) = fe(:);
    
    end
    
    phi = me\feall;
    
          sigmanp(1,:) = phi(:,1);
        
          sigmanp(2,:) = phi(:,2);

          sigmanp(3,:) = phi(:,3);

          sigmanp(4,:) = phi(:,4);  

          sigmanp(5,:) = phi(:,5);   

          sigmanp(6,:) = phi(:,6);

          epsilonnp(1,:) = phi(:,7);
        
          epsilonnp(2,:) = phi(:,8);
 
          epsilonnp(3,:) = phi(:,9);

          epsilonnp(4,:) = phi(:,10);  

          epsilonnp(5,:) = phi(:,11);   

          epsilonnp(6,:) = phi(:,12);



sigelnp = zeros(nnpe,1);

    for jnp = 1:1:nnpe
        sigelnp(jnp) = phi(jnp,maxcomp);
    end
    
 sigmaxnp = max(sigelnp);
 sigminnp = min(sigelnp);
 sigrangenpmax = sigmaxnp-sigminnp;
 sigrangenp(iele) = sigrangenpmax;
 
   if(sigrangenp(iele)  > 1.0*sigrangeqp(iele))
      damping_factors(iele) = damp_iter*0.01; 
   else
     dampflagel(iele) = 1; 
     dampflag = 1;
   end
 
 if(damp_iter > 20)
     dampflag = 1;
 end

%
end

for jnp = 1:1:nnpe
    j1 = newnp(iele,jnp);
    sigmanpave(:,j1) = sigmanpave(:,j1)+sigmanp(:,jnp);
    epsilonnpave(:,j1) = epsilonnpave(:,j1)+epsilonnp(:,jnp);
    numelepernp(j1) = numelepernp(j1)+1;
end
    


end

    for inp = 1:1:newnumnp
    sigmanpave(:,inp) = sigmanpave(:,inp)./numelepernp(inp);
    epsilonnpave(:,inp) = epsilonnpave(:,inp)./numelepernp(inp);
    end



    status = ' Finished stress smoothing routine ' 
    

   

%    next compute metrics related to stiffness (, ...)



    status = ' Starting calculation of stiffness-related metrics ' 

grain_inverse_rmatrix= zeros(6,6,ngrains);
embeddedstiffness = zeros(1,newnumnp);
re2sx = zeros(1,newnumnp);
SXStiffness = zeros(1,newnumnp);

for igrain = 1:ngrains
rotated_rmatrix = grain_rmatrix(:,:,igrain); 
inverse_rmatrix = inv(rotated_rmatrix);
grain_inverse_rmatrix(:,:,igrain) =inverse_rmatrix(:,:);
end

for inp = 1:1:newnumnp
    
    % look up grain number 
    
    igrain = grain4np(inp);
    
    inverse_rmatrix = grain_inverse_rmatrix(:,:,igrain);
    
%  component of the stiffness/compliance depends on loading mode
%    --  currently SX stiffness and RE2IS computed only for uniaxial extension modes
%    --  still thinking about what makes sense for these
%
   switch loadcode
       case 1
        embeddedstiffness(inp) = sigma_average(1)/epsilonnpave(1,inp); 
        SXStiffness(inp) = 1./inverse_rmatrix(1,1);
        re2sx(inp) = embeddedstiffness(inp)/SXStiffness(inp);
       case 2
        embeddedstiffness(inp) = sigma_average(2)/epsilonnpave(2,inp); 
        SXStiffness(inp) = 1./inverse_rmatrix(2,2);
        re2sx(inp) = embeddedstiffness(inp)/SXStiffness(inp);
       case 3   
        embeddedstiffness(inp) = sigma_average(3)/epsilonnpave(3,inp);
        SXStiffness(inp) = 1./inverse_rmatrix(3,3);
        re2sx(inp) = embeddedstiffness(inp)/SXStiffness(inp);
       case 4
        embeddedstiffness(inp) = (sigma_average(1)/epsilonnpave(1,inp)+sigma_average(2)/epsilonnpave(2,inp))/2.; 
%        re2sx(inp) = embeddedstiffness(inp)/rotated_rmatrix(1,1);
       case 5
        embeddedstiffness(inp) = (sigma_average(2)/epsilonnpave(2,inp)+sigma_average(3)/epsilonnpave(3,inp))/2.; 
%        re2sx(inp) = embeddedstiffness(inp)/rotated_rmatrix(1,1);
       case 6
        embeddedstiffness(inp) = (sigma_average(1)/epsilonnpave(1,inp)+sigma_average(3)/epsilonnpave(3,inp))/2.; 
%        re2sx(inp) = embeddedstiffness(inp)/rotated_rmatrix(1,1);
   end
 
end
%
%  now compute metrics that involve the SCYS
%
schmidfactors = zeros(1,newnumnp);
rse = zeros(1,newnumnp);
grain_schmidfactors = zeros(1,ngrains);
%
%        Compute Schmid tensor for each grain
%
 phase_schmid_tensors = PhaseSchmidTensor(nphases,crystal_type_all,covera_all,maxss);

 schmid_tensors = CalcSchmidTensors(rotations,phase_schmid_tensors,phaseofgrain,maxss);

%
%   now compute the strength-related metrics
%



    status = ' Starting calculation of (stiffness+strength)-related metrics ' 

for   inp =1:1:newnumnp
         
    
igrain = grain4np(inp);
schmid_tensor = schmid_tensors(:,:,igrain);
    
%pull elemental stress from global arrays
npsigma(:) = sigmanpave(:,inp);

%compute resolved shear stresses
rss = schmid_tensor*npsigma';

%compute strength-to-stiffness
elrse = max(rss./sss);
elrse = Eeff/elrse;

%compute schmid factor based on nominal stress
schmidfactor = max(abs(schmid_tensor*nomstress./sssr));

%save elemental quantities in global arrays

rse(inp) = elrse;
schmidfactors(inp) = schmidfactor;

%the next line redefines the schmidfactor over and over again, but they are
%all the same
grain_schmidfactors(igrain) = schmidfactor;
 
 
end

% repeat the computation of metrics for the element-averaged stress and
% strain

embeddedstiffness_ave = zeros(1,numel);
re2sx_ave = zeros(1,numel);
SXStiffness_ave = zeros(1,numel);
rse_ave= zeros(1,numel);
schmidfactors_ave = zeros(1,numel);
elsigma= zeros(6,1);
elepsilon= zeros(6,1);

for   iele =1:1:numel
         
igrain = grains(iele);
    
%pull elemental stress and strain from global arrays
elsigma(:) = sigma_ave(:,iele);
elepsilon(:) = epsilon_ave(:,iele);
    
inverse_rmatrix = grain_inverse_rmatrix(:,:,igrain);
 
   switch loadcode
       case 1
        embeddedstiffness_ave(iele) = sigma_average(1)/elepsilon(1); 
        SXStiffness_ave(iele) = 1./inverse_rmatrix(1,1);
        re2sx_ave(iele) = embeddedstiffness_ave(iele)/SXStiffness_ave(iele);
       case 2
        embeddedstiffness_ave(iele) = sigma_average(2)/elepsilon(2); 
        SXStiffness_ave(iele) = 1./inverse_rmatrix(2,2);
        re2sx_ave(iele) = embeddedstiffness_ave(iele)/SXStiffness_ave(iele);
       case 3   
        embeddedstiffness_ave(iele) = sigma_average(3)/elepsilon(3);
        SXStiffness_ave(iele) = 1./inverse_rmatrix(3,3);
        re2sx_ave(iele) = embeddedstiffness_ave(iele)/SXStiffness_ave(iele);
       case 4
        embeddedstiffness_ave(iele) = (sigma_average(1)/elepsilon(1)+sigma_average(2)/elepsilon(2))/2.; 
%        re2sx_ave(iele) = embeddedstiffness(inp)/rotated_rmatrix(1,1);
       case 5
        embeddedstiffness_ave(iele) = (sigma_average(2)/elepsilon(2)+sigma_average(3)/elepsilon(3))/2.; 
%        re2sx_ave(iee) = embeddedstiffness(inp)/rotated_rmatrix(1,1);
       case 6
        embeddedstiffness_ave(iele) = (sigma_average(1)/elepsilon(1)+sigma_average(3)/elepsilon(3))/2.; 
%        re2sx_ave(iele) = embeddedstiffness(inp)/rotated_rmatrix(1,1);
   end

schmid_tensor = schmid_tensors(:,:,igrain);

%compute resolved shear stresses
rss = schmid_tensor*elsigma;

%compute strength-to-stiffness
elrse = max(rss./sss);
elrse = Eeff/elrse;

%compute schmid factor based on nominal stress
schmidfactor = max(abs(schmid_tensor*nomstress./sssr));


%save elemental quantities in global arrays
rse_ave(iele) = elrse;
schmidfactors_ave(iele) = schmidfactor;

end

status = ' Finished caluclation of metrics  --  Starting Output to VTK file'

% map displacements to new mesh

ugbg = zeros(3*newnumnp,1);
for inp = 1:1:numnp
    ugbg(3*inp) = u(3*inp);
    ugbg(3*inp-1) = u(3*inp-1);
    ugbg(3*inp-2) = u(3*inp-2);
end
for inp = numnp+1:1:newnumnp
    inpold = old2new(inp); 
    ugbg(3*inp) = u(3*inpold);
    ugbg(3*inp-1) = u(3*inpold-1);
    ugbg(3*inp-2) = u(3*inpold-2);
end


solutionfilegbg    = Data_Structure_Solution(ugbg,sigma_ave,epsilon_ave,sigmanpave,epsilonnpave,rse,schmidfactors,embeddedstiffness,SXStiffness,re2sx,rse_ave,schmidfactors_ave,embeddedstiffness_ave,SXStiffness_ave,re2sx_ave);
graindatafile   = Data_Structure_GrainData(grain_volumes,rotations,grain_schmidfactors);

basenamegbg = input('Base name of grain-by-grain output VTK file: ', 's');

SUF = 'vtk';

fileVTKgbg = fopen(sprintf('%s.%s', [basenamegbg], SUF), 'w');

Export_2_VTK(fileVTKgbg, crystal_type, microstructurefile, meshfilegbg, solutionfilegbg,graindatafile,grain4np)

status = ' Finished writing output to VTK file'


ReturnStatus = 'Simulation Completed';
