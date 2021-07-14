function microstructurefile = MicrostructureStructure(phases,grains,orientations,rotations)
% MicrostructureStructure - Create microstructure structure from grain phase and orientation data.
%   
%   USAGE:
%
%   microstructurefile = MicrostructureStructure()
%
%   INPUT:
%
  %  the phase, grain, orientation and rotation arrays set up in ReadNeperMshFile.m

%   OUTPUT:
  %  
%    the matlab structure used in MechMet and Mechmonics
%        
%

microstructurefile = struct('phases',phases,'grains',grains,'orientations',orientations,'rotations',rotations)
