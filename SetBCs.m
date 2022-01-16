function [bccode,bcvalue,nomstress] = SetBCs(m,loadcode,nomstrain,neperversionnumber,varargin)
% SetBCs - Write boundary conditions based on loadcode parameter
%
%   INPUT:
%  
%  the code is defined as follows:
%   loadcode = 1  --  x-tension
%   loadcode = 2  --  y-tension
%   loadcode = 3  --  z-tension
%   loadcode = 4  --  x-y biaxial tension
%   loadcode = 5  --  y-z biaxial tension
%   loadcode = 6  --  z-x biaxial tension    
%
%     Boundary conditions are specified by a CODE and VALUE for each
%     of the three coordinate directions.  A BC specification is a
%     3 x 2 matrix.  The first column contains the CODES, nonzero for
%     an essential BC, zero for no essential BC.  The second column
%     gives the value of the applied displacement. 
%
%   OUTPUT:
%
%     bccode and bcvalue arrays with codes for applying bc's
%
%
   xcode =0;
   ycode =0;
   zcode =0;
   xval  =0;
   yval  =0;
   zval  =0;
   if(loadcode==1)
       xcode =1;
       xval = nomstrain(1);
       nomstress = [1;0;0;0;0;0];
   end
   if(loadcode==2)
       ycode =1;
       yval = nomstrain(2);
       nomstress = [0;1;0;0;0;0];
   end
   if(loadcode==3)
       zcode =1;
       zval = nomstrain(3);
       nomstress = [0;0;1;0;0;0];
   end
   if(loadcode==4)
       xcode =1;
       xval = nomstrain(1);     
       ycode =1;
       yval = nomstrain(2);
       nomstress = [1;1;0;0;0;0];
   end
   if(loadcode==5)
       ycode =1;
       yval = nomstrain(2);
       zcode =1;
       zval = nomstrain(3);
       nomstress = [0;1;1;0;0;0];
   end
   if(loadcode==6)
       zcode =1;
       zval = nomstrain(3);
       xcode =1;
       xval = nomstrain(1); 
       
       nomstress = [1;0;1;0;0;0];
   end

  optcell = {...
   'BCXMin',             [1 0; 0 0; 0 0;], ...
   'BCXMax',             [xcode xval; 0 0; 0 0;], ...
   'BCYMin',             [0 0; 1 0; 0 0;], ...
   'BCYMax',             [0 0; ycode yval; 0 0;], ...
   'BCZMin',             [0 0; 0 0; 1 0;], ...
   'BCZMax',             [0 0; 0 0; zcode zval;]  ...
      };

%
options = OptArgs(optcell, varargin);
%
%-------------------- *
%
if nargin < 3
  error(sprintf('Not enough args:\nUsage:  %s', Usage))
end
%
CODE = 1; VAL = 2;
%
%  Apply BC's.
%
%bcfile = sprintf('%s.bcs', fname);
%fprintf('boundary conditions file:\n         %s\n', bcfile);
%
surfopts = {'minx', 1, 'miny', 1, 'minz', 1, 'maxx', 1, 'maxy', 1, 'maxz', 1};

if(neperversionnumber==3)
surfs = RectMeshSurfaces(m, surfopts{:});
else
surfs = NeperMeshSurfaces(m, surfopts{:});
end

surfspecs = struct('name',   {'minx',   'miny',   'minz',   'maxx',   'maxy',   'maxz'},  ...
		   'bcname', {'BCXMin', 'BCYMin', 'BCZMin', 'BCXMax', 'BCYMax', 'BCZMax'} ...
		   );  % to be continued ...2006-02-28  16:18
%
nnps = size(m.crd, 2);
bc_codes = zeros(3, nnps);
bc_vals  = zeros(3, nnps);
%
for spec=surfspecs
  numpts = length(surfs.(spec.name));
  i = repmat([1;2;3], [1 numpts]);
  j = repmat(surfs.(spec.name), [3 1]);
  v = repmat(options.(spec.bcname)(:, CODE), [1 numpts]);
  bc_codes = bc_codes + sparse(i(:), j(:), v(:), 3, nnps);
  v = repmat(options.(spec.bcname)(:, VAL), [1 numpts]);
  bc_vals  = bc_vals  + sparse(i(:), j(:), v(:), 3, nnps);
end
%
bc_codes(bc_codes ~=0) = 1;
bcnz = find(sum(bc_codes));
nbc  = length(bcnz);
%
 vfmt = '%20.12e %20.12e %20.12e\n';
for node=bcnz
  bc = bc_codes(:, node)';
  bv = bc_vals(:, node)';
  sc = [(node - 1), bc];
  sv = [(node - 1), bv];
%  fprintf(f, s);
  if (node==bcnz(1))
   sctotal = sc;
   svtotal = sv;
  else
   sctotal = cat(1,sctotal,sc);
   svtotal = cat(1,svtotal,sv);
  end
end
%
% fprintf(f, '0\n');
% fclose(f);
bccode  = sctotal;
bcvalue = svtotal;

%
%--------------------*--------------------------------------------------
%
function s = TrueFalse(code)
% TRUEFALSE -
%
FT   = 'FT';
code = code + 1;
sp = '   ';
%
s = [sp, FT(code(1)), sp, FT(code(2)), sp, FT(code(3)), sp];
%
