function [L, bX, bY] = NuclearSplitting(NuclearMask, TauMinArea,...
                                        TauMaxArea)
%Performs splitting of clumped nuclei using watershed and h-minima
%transform to generate label mask.
%
%inputs:
%NuclearMask - (M x N logical) Mask of nuclear pixels.
%TauMinArea - (scalar) minimum area in pixels for a single object to 
%               include in mask. Default value 30.
%TauMaxArea - (scalar) maximum area in pixels for a single object to 
%               include in mask. Default value 1000.
%
%outputs:
%L - (M x N float) Nuclei label image.
%bX - (K-length cell double) Cell array of doubles describing horizontal 
%      object boundary coordinates.
%bY - (K-length cell double) Cell array of doubles describing vertical
%      object boundary coordinates.
%
%notes:
%Kong J, Cooper LAD, et al "Machine-based morphologic analysis of 
%glioblastoma using whole-slide pathology images uncovers clinically 
%relevant molecular correlates," PLoS One. 2013 Nov 13;8(11):e81049. 
%doi: 10.1371/journal.pone.0081049. eCollection 2013.
%
%Authors: Jun Kong and Lee Cooper, Emory University.


%Parse inputs and set default values
switch nargin
    case 1
        TauMinArea = 30;
        TauMaxArea = 1000;
    case 2
        TauMaxArea = 1000;
end

%Generate distance, h-minima transforms of nuclear mask
Distance = -bwdist(~NuclearMask);
Distance(~NuclearMask) = -Inf;
HMinima = imhmin(Distance, 1);
NuclearMask(watershed(HMinima) == 0) = 0;
clear Distance HMinima;

%Zero non-minima and filter on nuclear area
L = bwlabel(NuclearMask, 4);
Area = regionprops(L, 'Area');
Area = cat(1, Area.Area);
Keep = find(Area > TauMinArea & Area < TauMaxArea);
AreaMask = ismember(L, Keep);
AreaMask = imfill(AreaMask, 'holes');

%Generate final label image and boundaries
L = bwlabel(AreaMask,4);
Bounds = bwboundaries(AreaMask, 4, 'noholes');
bX = cellfun(@(x)x(:,2), Bounds, 'UniformOutput', false);
bY = cellfun(@(x)x(:,1), Bounds, 'UniformOutput', false);
