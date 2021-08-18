function surfs = NeperMeshSurfaces(m, varargin)
% NeperMeshSurfaces - Surfaces of a Neper generated mesh.
%
%   USAGE:
%
%   surfs = NeperMeshSurfaces(m, <options>)
%
%   INPUT:
%
%  m is a MeshStructure 
%    It delineates a six-sided region.
%
%   <options> is a sequence of parameter, value pairs
%             Available options are listed below with default values
%             shown in brackets.
%
%    'minx'      {0} | 1
%    'miny'      {0} | 1
%    'minz'      {0} | 1
%    'maxx'      {0} | 1
%    'maxy'      {0} | 1
%    'maxz'      {0} | 1
%
%   OUTPUT:
%
%   surfs is a structure
%         It has components 'minx', etc., for each surface turned
%         on by the options.
%
%   NOTES:
%
%   * Each componenet of surfs is currently just a list of node numbers
%   
MyName = mfilename;
%
%-------------------- Defaults and Keyword Arguments
%
optcell = {...
    'minx',  0, ...
    'miny',  0, ...
    'minz',  0, ...
    'maxx',  0, ...
    'maxy',  0, ...
    'maxz',  0
       };
%
opts = OptArgs(optcell, varargin);
%
%-------------------- *
%
%  Need to find surface elements on boundary and number of 
%  solid element containing it.
%
%

ELSURFS = m.etype.surfs;
%
sz = size(ELSURFS); nvels = size(m.con, 2); nsels = nvels*sz(2);
%
SURFELS = reshape(m.con(ELSURFS, :), [sz(1), nsels]);
VOLELS  = reshape(repmat(1:nvels, [sz(2), 1]), [1, sz(2)*nvels]);
%
if opts.minx
    nse = size(m.nps_1,1);
    np = m.nps_1(:,[2:7]);
    surfs.minx = unique(np)';
end
%
if opts.miny
    nse = size(m.nps_3,1);
    np = m.nps_3(:,[2:7]);
    surfs.miny = unique(np)';
end
%
if opts.minz
    nse = size(m.nps_5,1);
    np = m.nps_5(:,[2:7]);
    surfs.minz = unique(np)';
end
%
if opts.maxx
    nse = size(m.nps_2,1);
    np = m.nps_2(:,[2:7]);
    surfs.maxx = unique(np)';
end
%
if opts.maxy
    nse = size(m.nps_4,1);
    np = m.nps_4(:,[2:7]);
    surfs.maxy = unique(np)';
end
%
if opts.maxz
    nse = size(m.nps_6,1);
    np = m.nps_6(:,[2:7]);
    surfs.maxz = unique(np)';
end
%
