function [Features, Names, cX, cY] = FeatureExtraction(L, I, K, FSDBins,...
                                                        Delta, M)
%Extract shape, texture, gradient and intensity features from segmented 
%objects.
%
%inputs:
%L - (T x T float) Nuclei label image.
%I - (T x T x 3 uint 8) Color image.
%K - (scalar) Points for boundary resampling in calculating fourier 
%    descriptors. Default value 128.
%FSDBins - (scalar) number of frequency bins for calculating fourier shape
%           descriptors. Default value 6.
%Delta - (scalar) Dilation factor for expanding nuclei to capture
%        surrounding cytoplasm.
%M - (3 x 3 float) Color deconvolution matrix for deconvolving hematoxlyin
%    and eosin signals. Default value [0.650 0.072 0; 0.704 0.990 0; 0 0 0]
%
%outputs:
%Features - (N x 48 float) Matrix of features where each row is a nucleus
%           and each column is a feature.
%Names - (48-length cell string) Cell array of strings describing features.
%cX - (N x 1 float) Vector of horizontal nuclear centroid coordinates.
%cY - (N x 1 float) Vector of vertical nuclear centroid coordinates.
%
%Notes:
%Kong J, Cooper LAD, et al "Machine-based morphologic analysis of 
%glioblastoma using whole-slide pathology images uncovers clinically 
%relevant molecular correlates," PLoS One. 2013 Nov 13;8(11):e81049. 
%doi: 10.1371/journal.pone.0081049. eCollection 2013.
%
%Authors: Lee Cooper and Jun Kong, Emory University.


%Parse inputs and set default values
switch nargin
    case 2
        K = 128;
        FSDBins = 6;
        Delta = 8;
        M = [0.650 0.072 0; 0.704 0.990 0; 0 0 0];
    case 3
        FSDBins = 6;
        Delta = 8;
        M = [0.650 0.072 0; 0.704 0.990 0; 0 0 0];
    case 4
        Delta = 8;
        M = [0.650 0.072 0; 0.704 0.990 0; 0 0 0];
    case 5
        M = [0.650 0.072 0; 0.704 0.990 0; 0 0 0];
end

%Get number of objects in label mask
N = max(L(:));

%Built-in shape features                    
statsI = regionprops(L, 'Area','Perimeter','Eccentricity',...
                        'MajorAxisLength','MinorAxisLength','Extent',...
                        'Solidity','PixelIdxList','Centroid',...
                        'BoundingBox');
fArea = cat(1,statsI.Area);
fPerimeter = cat(1,statsI.Perimeter);
fEccentricity = cat(1,statsI.Eccentricity);
fCircularity = 4*pi * fArea./ (fPerimeter.^2);
fMajorAxisLength = cat(1,statsI.MajorAxisLength);
fMinorAxisLength = cat(1,statsI.MinorAxisLength);
fExtent = cat(1,statsI.Extent);
fSolidity = cat(1,statsI.Solidity);
fMorph = [fArea, fPerimeter, fEccentricity, fCircularity,...
            fMajorAxisLength, fMinorAxisLength, fExtent, fSolidity];
MorphNames = {'Area', 'Perimeter', 'Eccentricity', 'Circularity',...
                'MajorAxisLength', 'MinorAxisLength', 'Extent',...
                'Solidity'};

%Unpack nuclear centroids
Centroids = cat(1, statsI.Centroid);
cX = Centroids(:, 1);
cY = Centroids(:, 2);
            
%Generate object pixel lists for nuclear and cytoplasmic regions
Nuclei = cell(1, N);
Cytoplasms = cell(1, N);
Bounds = cell(1, N);
disk = strel('disk', Delta, 0); %create round structuring element
for i = 1:N
    Nuclei{i} = statsI(i).PixelIdxList;
    bounds = GetBounds(statsI(i).BoundingBox, Delta, size(L,1), size(L,2));
    Nucleus = L(bounds(3):bounds(4), bounds(1):bounds(2)) == i;
    Trace = bwboundaries(Nucleus, 8, 'noholes');
    Bounds{i} = Trace{1};
    Mask = L(bounds(3):bounds(4), bounds(1):bounds(2)) > 0;
    cytoplasm = xor(Mask, imdilate(Nucleus, disk));
    Cytoplasms{i} = PixIndex(cytoplasm, bounds, size(L,1), size(L,2));
end
            
%Calculate Fourier shape descriptors
Interval = round(Log2Spaced(0, log2(K)-1, FSDBins+1));
FSDNames = cellfun(@(x,y) [x num2str(y)], repmat({'FSD'}, [1,FSDBins]),...
                    num2cell(1:FSDBins), 'UniformOutput', false);
FSDGroup = zeros(N, FSDBins);
for i = 1:N
    FSDGroup(i,:) = FourierShapeDescriptors(Bounds{i}(:,1),...
        Bounds{i}(:,2), K, Interval);
end

%Deconvolve color image to calculate nuclear, cytoplasmic texture features
Deconvolved = ColorDeconvolution(I, M, [true true false]);
Hematoxylin = Deconvolved(:,:,1);
Hematoxylin(Hematoxylin > 255) = 255;
Eosin = Deconvolved(:,:,2);
Eosin(Eosin > 255) = 255;
clear Deconvolved;

%Convert deconvolved images to double format for feature extraction
Hematoxylin = double(Hematoxylin);
Eosin = double(Eosin);

%Hematoxlyin features calculation in nuclear regions
[HematoxylinIntensityGroup, IntensityNames] = ...
    IntensityFeatureGroup(Hematoxylin, Nuclei);
[HematoxylinTextureGroup, TextureNames] = ...
    TextureFeatureGroup(Hematoxylin, Nuclei);
[HematoxylinGradientGroup, GradientNames] = ...
    GradientFeatureGroup(Hematoxylin, Nuclei);

%Eosin feature calculation in cytoplasm regions
EosinIntensityGroup = IntensityFeatureGroup(Eosin, Cytoplasms);
EosinTextureGroup = TextureFeatureGroup(Eosin, Cytoplasms);
EosinGradientGroup = GradientFeatureGroup(Eosin, Cytoplasms);

%Combine feature sets
fDeconvolved = [HematoxylinIntensityGroup HematoxylinTextureGroup...
    HematoxylinGradientGroup EosinIntensityGroup EosinTextureGroup...
    EosinGradientGroup];

%concatenate features
Features = [fMorph FSDGroup fDeconvolved];
   
%Generate feature names output
Names = [IntensityNames TextureNames GradientNames];
NuclearNames = cellfun(@(x)strcat('Hematoxlyin', x), Names,...
    'UniformOutput', false);
CytoplasmNames = cellfun(@(x)strcat('Cytoplasm', x), Names, 'UniformOutput', false);
Names = [MorphNames FSDNames NuclearNames CytoplasmNames];

end


function [Intensity, Names] = IntensityFeatureGroup(I, ObjectPixelList)
Intensity = zeros(length(ObjectPixelList), 4);
for i = 1:length(ObjectPixelList)
    pixOfInterest = I(ObjectPixelList{i});
    Intensity(i,1) = double(mean(pixOfInterest));
    Intensity(i,2) = Intensity(i,1) - double(median(pixOfInterest));
    Intensity(i,3) = max(pixOfInterest);
    Intensity(i,4) = min(pixOfInterest);
    Intensity(i,5) = std(double(pixOfInterest));
end
Names = {'MeanIntensity', 'MeanMedianDifferenceIntensity',...
    'MaxIntensity', 'MinIntensity', 'StdIntensity'};
end

function [Features, Names]= TextureFeatureGroup(I, ObjectPixelList)
Features = zeros(length(ObjectPixelList), 4);
for i = 1:length(ObjectPixelList)
    pixOfInterest = I(ObjectPixelList{i});
    [counts] = imhist(uint8(pixOfInterest));
    prob = counts/sum(counts);
    Features(i,1) = entropy(uint8(pixOfInterest));
    Features(i,2) = sum(prob.^2);
    Features(i,3) = skewness( double(pixOfInterest) );
    Features(i,4) = kurtosis( double(pixOfInterest) );
end
Names = {'Entropy', 'Energy', 'Skewness', 'Kurtosis'};
end

function [Features, Names] = GradientFeatureGroup(I, ObjectPixelList)
[Gx, Gy] = gradient(double(I));
diffG = sqrt(Gx.*Gx+Gy.*Gy);
BW_canny = edge(I,'canny');

Features = zeros(length(ObjectPixelList), 8);
for i = 1:length(ObjectPixelList)
    pixOfInterest = diffG(ObjectPixelList{i});
    fMeanGradMag = mean(pixOfInterest);
    fStdGradMag = std(pixOfInterest);
    [counts, ~] = imhist(uint8(pixOfInterest));
    prob = counts/sum(counts);
    
    fEntropyGradMag = entropy(uint8(pixOfInterest));
    fEnergyGradMag = sum(prob.^2);
    fSkewnessGradMag = skewness( double(pixOfInterest) );
    fKurtosisGradMag = kurtosis( double(pixOfInterest) );
    
    bw_canny = BW_canny(ObjectPixelList{i});
    fSumCanny = sum(bw_canny(:));
    
    fMeanCanny = fSumCanny / length(pixOfInterest);
    
    Features(i,:) = [fMeanGradMag, fStdGradMag, fEntropyGradMag,...
        fEnergyGradMag, fSkewnessGradMag,fKurtosisGradMag, fSumCanny,...
        fMeanCanny];
end
Names = {'MeanGradMag', 'StdGradMag', 'EntropyGradMag', 'EnergyGradMag',...
    'SkewnessGradMag', 'KurtosisGradMag', 'SumCanny', 'MeanCanny'};
end

function bounds = GetBounds(bbox, delta, M, N)
%get bounds of object in global label image
bounds(1) = max(1,floor(bbox(1) - delta));
bounds(2) = min(N, ceil(bbox(1) + bbox(3) + delta));
bounds(3) = max(1,floor(bbox(2) - delta));
bounds(4) = min(M, ceil(bbox(2) + bbox(4) + delta));
end

function idx = PixIndex(Binary, bounds, M, N)
%get global linear indices of object extracted from tile
[i, j] = find(Binary);
i = i + bounds(3) - 1;
j = j + bounds(1) - 1;
idx = sub2ind([M N], i, j);
end
