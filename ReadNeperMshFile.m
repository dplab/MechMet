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

NodeFlag = '$Nodes';
ElementFlag = '$Elements';
OrientationFlag = '$ElsetOrientations';
GroupFlag = '$Groups';
FasetsFlag = '$Fasets';

x0 = 'x0';
x1 = 'x1';
y0 = 'y0';
y1 = 'y1';
z0 = 'z0';
z1 = 'z1';

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
       numnp = str2num(fgetl(fid)) ;
       crd= zeros(3,numnp);
       tmp_coords = fscanf(fid, '%d %f %f %f', 4*numnp);
       tmp_coords = reshape(tmp_coords, [4, numnp]);
       crd = tmp_coords(2:4,:);
       clear tmp_coords;
    end
 %
    if(tstring== ElementFlag) 
       tmp_numel = str2num(fgetl(fid)) ;
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
       numori = tmp_ori1;
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
    
       if(tstring== FasetsFlag)  
           
       nblocks = str2num(fgetl(fid));
       
       for nbl = 1:1:nblocks
           
       minmax_flag = convertCharsToStrings(fgetl(fid));   
      
       if(minmax_flag == x0)           
       nse_1 = str2num(fgetl(fid)) ;
%       np_tmp1 = zeros(nse_1,7);
       nps_1 = fscanf(fid, '%d %d %d %d %d %d %d', 7*nse_1);
       nps_1 = reshape(nps_1, [ 7 , nse_1])';
       dum = fgetl(fid);
       end
       
       if(minmax_flag == x1)
       nse_2 = str2num(fgetl(fid)) ;
%       np_tmp2 = zeros(nse_2,7);
       nps_2 = fscanf(fid, '%d %d %d %d %d %d %d', 7*nse_2);
       nps_2 = reshape(nps_2, [7,  nse_2])';
       dum = fgetl(fid);
       end
       
       if(minmax_flag == y0)
       nse_3 = str2num(fgetl(fid)) ;
%       np_tmp3 = zeros(nse_3,7);
       nps_3 = fscanf(fid, '%d %d %d %d %d %d %d', 7*nse_3);
       nps_3 = reshape(nps_3, [7, nse_3])';
       dum = fgetl(fid);
       end
       
       if(minmax_flag == y1)
       nse_4 = str2num(fgetl(fid)) ;
%       np_tmp4 = zeros(nse_4,7);
       nps_4 = fscanf(fid, '%d %d %d %d %d %d %d', 7*nse_4);
       nps_4 = reshape(nps_4, [7, nse_4])';
       dum = fgetl(fid);
       end
             
       if(minmax_flag == z0)        
       nse_5 = str2num(fgetl(fid)) ;
%       np_tmp5 = zeros(nse_5,7);
       nps_5 = fscanf(fid, '%d %d %d %d %d %d %d', 7*nse_5);
       nps_5 = reshape(nps_5, [7, nse_5])';
       dum = fgetl(fid);
       end
       
       if(minmax_flag == z1)              
       nse_6 = str2num(fgetl(fid)) ;
%       np_tmp6 = zeros(nse_6,7);
       nps_6 = fscanf(fid, '%d %d %d %d %d %d %d', 7*nse_6);
       nps_6 = reshape(nps_6, [7, nse_6])';
       dum = fgetl(fid);
       end
       
       end
       
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

meshfile.nps_1 = nps_1;
meshfile.nps_2 = nps_2;
meshfile.nps_3 = nps_3;
meshfile.nps_4 = nps_4;
meshfile.nps_5 = nps_5;
meshfile.nps_6 = nps_6;

microstructurefile = MicrostructureStructure(grains(2,:),grains(1,:),orientations,rotations);
microstructurefile.name = [fname '_micro'];
%
fclose(fid);

