function NuclearMask = NuclearMask(I, TauMinArea, TauMaxArea,...
                                        TauRedStrict, TauRedRelaxed,...
                                        TauMorph, TauNuclei)
%Nuclei foreground segmentation using morphological operations. Bi-level
%thresholding is used to increase detection sensitivity while minimizing
%noise. Returns a binary mask with objects in range Area \in [TauMinArea, 
%TauMaxArea].
%
%inputs:
%I - (T x T x 3 uint8) RGB color image.
%TauMinArea - (scalar) minimum area in pixels for a single object to 
%               include in mask. Default value 30.
%TauMaxArea - (scalar) maximum area in pixels for a single object to 
%               include in mask. Default value 1000.
%TauRedStrict - (scalar) strict threshold on red/green ratio used to 
%               identify candidate red blood cell pixels. Default value 5.
%TauRedRelaxed - (scalar) relaxed threshold on red/green ratio used to 
%               identify candidate red blood cell pixels. Default value 4.
%TauMorph - (scalar) strict threshold used to identify foreground in 
%           morphological reconstruction.
%TauNuclei - (scalar) relaxed threshold used to identify foreground in 
%           morphological reconstruction.
%
%outputs:
%NuclearMask - (T x T logical) image where nuclear pixels have value true.
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
        TauRedStrict = 5;
        TauRedRelaxed = 4;
        TauMorph = 80;
        TauNuclei = 45;
    case 2
        TauMaxArea = 1000;
        TauRedStrict = 5;
        TauRedRelaxed = 4;
        TauMorph = 80;
        TauNuclei = 45;
    case 3
        TauRedStrict = 5;
        TauRedRelaxed = 4;
        TauMorph = 80;
        TauNuclei = 45;
    case 4
        TauRedRelaxed = 4;
        TauMorph = 80;
        TauNuclei = 45;
    case 5
        TauMorph = 80;
        TauNuclei = 45;
    case 6
        TauNuclei = 45;
end

%Mask red blood cells using color information
RedRatio = single(I(:,:,1))./(single(I(:,:,2))+eps);
[Y, X] = find(RedRatio > TauRedStrict);
if ~isempty(X)
    RBCMask = bwselect(RedRatio > TauRedRelaxed, X, Y, 8) & ...
        (double(I(:,:,1))./(double(I(:,:,3))+eps)>1);
else
    RBCMask = false(size(RedRatio));
end

%Remove background noise with morphological reconstruction
Negated = 255 - I(:,:,1);
Reconstructed = imreconstruct(imopen(Negated, strel('disk',10)), Negated);
Difference = Negated-Reconstructed;

%Filter on object area
Foreground = imfill(Difference > TauMorph, 'holes'); %strict thresholding
L = bwlabel(Foreground, 8);
Area = regionprops(L, 'Area');
Area = cat(1, Area.Area);
Keep = find(Area > TauMinArea & Area < TauMaxArea);
KeepMask = ismember(L, Keep);
[Y, X] = find(KeepMask);

%Process further is viable objects present
if ~isempty(X)

    %Threshold reconstruction (relaxed)
    NuclearMask = Difference > TauNuclei;
    
    %Select objects passing area threshold
    NuclearMask = bwselect(NuclearMask, X, Y, 8) & ~RBCMask;
    
    %Fill holes
    NuclearMask = imfill(NuclearMask, 'holes');
    
    %Open to break thin connections
    NuclearMask = imopen(NuclearMask, strel('disk',1));
    
    %Remove and debris produced by opening
    NuclearMask = imdilate(bwareaopen(NuclearMask,30),strel('disk',1));

else %Return blank image

    NuclearMask = false(size(RedRatio));

end
