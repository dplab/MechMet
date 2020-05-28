function m = ReadNeperMesh(fname)
% ReadNeperMesh - read neper mesh from a file
%
%   STATUS: in development: usable
%
%   USAGE:
%
%   m = ReadNeperMesh(fname, <options>)
%
%   INPUT:
%
%   fname is a string for the full name of the mesh files to be read;

%
%   OUTPUT:
%
%   m is a structure (see MeshStructure)
%
%   NOTES:
%
%   * Currently only 10-node tets are available
%   * Currently, connectivity is returned in fepx order
%
meshfile = [fname,'.mesh'];
%
fid = fopen(meshfile, 'r');
%
%  Header line
%
parms = fscanf(fid, '%d', 3);
numel = parms(1);
numnp = parms(2);
perel = parms(3) + 1;
%
%  Read connectivity.
%
con = fscanf(fid, '%d', (perel+1)*numel);
con = reshape(con, [(perel+1), numel]);
con = con(2:end, :) + 1;
%
%  Read coordinates.
%
tmp = fscanf(fid, '%d %f %f %f', 4*numnp);
tmp = reshape(tmp, [4, numnp]);

%

crd = tmp(2:4, :);
%
% *** To Be Done: read surface information
%
fclose(fid);
%
m = MeshStructure(crd, con, [], 'ElementType', 'tets:10');
m.name = fname;
%


