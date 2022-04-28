function [op,lp,ln] = layout_sdl(filename)
%LAYOUT_SDL Extract and plot optode and landmark points from a layout file
%
% [op, lp, ln] = LUMOFILE.UTIL.LAYOUT_SDL(filename) 
%
% LUMOFILE.UTIL.LAYOUT_SDL loads a LUMO layout file, extracts the 3D locations of each
% optode, landmark, and the landmark names. 
%
%   Parameters
%
%   fn: input filename
%
%   Returns
%
%   op: optode points, a <no. optodes x 3> matrix of the 3D locations of the optodes
%   lp: landmark points, a <no. landmarks x 3> matrix of the 3D locations of the landmarks
%   ln: landmark names, a cell array of landmark names for each entry in lp
%
%
%   (C) Gowerlabs Ltd., 2022
%

layout = lumofile.read_layout(filename);

optodes = [layout.docks.optodes];
no = length(optodes);
op = zeros(no, 3);

for i = 1:no
  op(i,:) = [optodes(i).coords_3d.x optodes(i).coords_3d.y optodes(i).coords_3d.z];
end

landmarks = [layout.landmarks];
nl= length(landmarks);
lp = zeros(nl, 3);

for i = 1:nl
  lp(i,:) = [landmarks(i).coords_3d.x landmarks(i).coords_3d.y landmarks(i).coords_3d.z];
end

ln = cell(nl, 1);
for i = 1:nl
  ln{i} = landmarks(i).name;
end
  
figure()
scatter3(op(:,1), op(:,2), op(:,3), 'o');
hold on;
scatter3(lp(:,1), lp(:,2), lp(:,3), 'x')
text(lp(:,1), lp(:,2), lp(:,3), ln, 'FontSize', 12);
hold off;
axis equal;
xlabel('x');
xlabel('y');
xlabel('z');
end

