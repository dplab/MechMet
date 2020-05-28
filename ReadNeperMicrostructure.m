function [g, p, k] = ReadNeperMicrostructure(micro)
% ReadNeperMicrostructure - Read microstructure information written by
% Neper
%
%
%   [g, p, k] = ReadNeperMicrostructure(micro)
%
%   INPUT:
%
%   micro is a string
%         the basename of the grains/orientations files
%
%   OUTPUT:
%
%   g is 1 X n (int)
%     the grain ID of each element
%   p is 1 X n (int)
%     the phase ID of each element
%   k is 3 X ng (float)
%     the kocks angles in degrees of each grain
%
%
% ===== Read grains and phases
%
grainfile = sprintf('%s.grain', micro);

f = fopen(grainfile, 'r');
%
tmp = fscanf(f,'%d %d', 2);
numel = tmp(1);
ngrn  = tmp(2);
%
%  * grain and phase data
%
tmp = fscanf(f,'%d %d', [2 numel]);
g = tmp(1, :);
p = tmp(2, :);
fclose(f);
%
% ===== Read orientations
%
orifile = sprintf('%s.kocks', micro);

f = fopen(orifile, 'r');
tmp = fgetl(f); % header line
ng = fscanf(f, '%d', 1);
tmp = fscanf(f, '%f', [4 ng]);
k = tmp(1:3, :);
fclose(f);
