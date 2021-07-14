function [meshfile,microstructurefile] = ReadNeperMshFile(fname)
% ReadNeperMshFile - read neper mesh, grains phases and orientations from a V4 *.msh file
%
%   STATUS: in development: usable
%
%   USAGE:
%
%   m = ReadNeperMshFile(fname, <options>)
%
%   INPUT:
%
%   fname is a string for the name of the msh files to be read (without the .msh extension);

%
%   OUTPUT:
%
%   meshfile and microstruturefile matlab structures (see MeshStructure)
%
%   NOTES:
%
%   * Currently only 10-node tets are available
%   * Currently, connectivity is returned in fepx order
%

NodeFlag = '$Nodes'
ElementFlag = '$Elements'
OrientationFlag = '$ElsetOrientations'
GroupFlag = '$Groups'

element_count=0;
phase_count = 1;
meshfile = [fname,'.msh'];
%
fid = fopen(meshfile, 'r');

while ~feof(fid)
    tline = fgetl(fid);
    tstring = convertCharsToStrings(tline);
 %
    if(tstring== NodeFlag)
       numnp = str2num(fgetl(fid)) 
       crd= zeros(3,numnp);
       tmp_coords = fscanf(fid, '%d %f %f %f', 4*numnp);
       tmp_coords = reshape(tmp_coords, [4, numnp]);
       crd = tmp_coords(2:4,:);
       clear tmp_coords;
    end
 %
    if(tstring== ElementFlag) 
       tmp_numel = str2num(fgetl(fid)) 
       tmp_np= zeros(10,tmp_numel);
       tmp_grains = zeros(1,tmp_numel);
       for iele = 1:1:tmp_numel
       tmp_element = str2num(fgetl(fid));
       if(tmp_element(2)==11)
           element_count = element_count+1;
           tmp_np([1:10],iele) = tmp_element([7:16]);
           tmp_grains(iele) = tmp_element(4);
       end
       end
    end
 %    
      if(tstring== OrientationFlag)
          
       tmp_ori1 = fscanf(fid, '%d', 1);
       tmp_ori2 = fgetl(fid);
       numori = tmp_ori1
       orientations= zeros(3,numori);
       tmp_oris = fscanf(fid, '%d %f %f %f', 4*numori);
       tmp_oris = reshape(tmp_oris, [4, numori]);
      end
  %
      if(tstring== GroupFlag)
       phase_count =2;
       tmp_group = fgetl(fid);
       numgrain = fscanf(fid, '%d', 1);
       phases = zeros(2,numgrain);
       tmp_phase = fscanf(fid, '%d %d', 2*numgrain);
       tmp_phase = reshape(tmp_phase, [2, numgrain]);
    end
 %
end

numel = element_count;
con = zeros(numel,10);
grains=zeros(2,numel);
orientations=zeros(3,numori);
icount = 1;
for iele=1:1:tmp_numel
    if(tmp_np(1,iele)~=0)
        con(icount,[1 3 5 10 2 4 6 7 9 8]) = tmp_np([1 2 3 4 5 6 7 8  9 10],iele);
        grainnum = tmp_grains(iele);
        grains(1,icount) = grainnum;
        if(phase_count==1)
            grains(2,icount) = 1;
        else
            grains(2,icount) = tmp_phase(2,grainnum);
        end
        icount = icount+1;
    end  
end
for iori = 1:1:numori
 orientations(:,iori) = tmp_oris(2:4,iori);    
end
quaterions = QuatOfRod(orientations);
rotations  = RMatOfQuat(quaterions);

meshfile = MeshStructure(crd, con', [], 'ElementType', 'tets:10');
meshfile.name = [fname '_mesh'];

microstructurefile = MicrostructureStructure(grains(2,:),grains(1,:),orientations,rotations);
microstructurefile.name = [fname '_micro'];
%
fclose(fid);

