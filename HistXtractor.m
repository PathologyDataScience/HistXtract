function HistXtractor(InputFile, Objective, T, OutputFolder,...
                            GenerateText, GenerateImages,...
                            TargetMu, TargetSigma, W)
%Extracts nuclear morphometric features from whole-slide images. Openslide
%is used to tile the input slide and each slide is analyzed independently.
%This analysis generates a 48-feature vector describing the shape, texture 
%and color information for each nucleus, along with centroid and boundary
%information.

%inputs:
%InputFile - (string) Full path and filename to slide image with extension.
%Objective - (scalar) Desired magnification for analysis.
%T - (scalar) Tile size in pixels for tiled processing of whole slide
%    image.
%OutputFolder - (string) Path where outputs will be generated.
%GenerateText - (logical) True indicates that text outputs for database
%               consumption should be produced.
%GenerateImages - (logical) True indicates that image outputs showing
%                 segmentations will be generated.
%TargetMu - (3-length float) Mean values of target color normalization
%           image in LAB color space. Can be generated from a sample image
%           tile where color deconvolution performs well. Default value
%           [-0.615876637607728 -0.043269404801974 0.038189608999642].
%TargetSigma - (3-length float) Mean values of target color normalization
%              image in LAB color space. Can be generated from a sample 
%              image tile where color deconvolution performs well. Default 
%              value [-0.615876637607728 -0.043269404801974 
%              0.038189608999642].
%W - (2 x 4 float) Linear discriminant parameters for masking tissue
%    from background. Used for color normalization. Default value
%    [-0.154 0.035 0.549 -45.718; -0.057 -0.817 1.170 -49.887].
%
%outputs:
%*** This function generates .mat files containing nuclear features and
%boundary and centroid information. Additional input switches enable the
%generation of .txt formats of this information for database ingestion, and
%image outputs for visualizing segmentation.
%
%Notes:
%nested function inputs:
%*** Nested functions contain a large number of parameters. These are not
%exposed in this function for brevity.
%NuclearMask: TauMinArea, TauMaxArea, TauRedStrict, TauRedRelaxed,
%               TauMorph, TauNuclei
%NuclearSplitting: NuclearMask, TauMinArea, TauMaxArea
%FeatureExtraction: L, I, K, FSDBins, Delta, M
%
%Kong J, Cooper LAD, et al "Machine-based morphologic analysis of 
%glioblastoma using whole-slide pathology images uncovers clinically 
%relevant molecular correlates," PLoS One. 2013 Nov 13;8(11):e81049. 
%doi: 10.1371/journal.pone.0081049. eCollection 2013.
%
%Authors: Lee Cooper and Jun Kong, Emory University.


%Parse inputs and set default values
switch nargin
    case 4
        GenerateText = false;
        GenerateImages = false;
        TargetMu = [-0.615876637607728 -0.043269404801974 ...
            0.038189608999642].';
        TargetSigma = [0.256490877567354 0.053757088006721 ...
            0.011683678645071].';
        W = [-0.154 0.035 0.549 -45.718; -0.057 -0.817 1.170 -49.887];
    case 5
        GenerateImages = false;
        TargetMu = [-0.615876637607728 -0.043269404801974 ...
            0.038189608999642].';
        TargetSigma = [0.256490877567354 0.053757088006721 ...
            0.011683678645071].';
        W = [-0.154 0.035 0.549 -45.718; -0.057 -0.817 1.170 -49.887];        
    case 6
        TargetMu = [-0.615876637607728 -0.043269404801974 ...
            0.038189608999642].';
        TargetSigma = [0.256490877567354 0.053757088006721 ...
            0.011683678645071].';
        W = [-0.154 0.035 0.549 -45.718; -0.057 -0.817 1.170 -49.887];
    case 7
        TargetSigma = [0.256490877567354 0.053757088006721 ...
            0.011683678645071].';
        W = [-0.154 0.035 0.549 -45.718; -0.057 -0.817 1.170 -49.887];
    case 8
        W = [-0.154 0.035 0.549 -45.718; -0.057 -0.817 1.170 -49.887];
end

%Add paths to helper functions
Path = mfilename('fullpath');
Slash = strfind(Path, '/');
if ~isempty(Slash)
    Path = Path(1:Slash(end));
end
% addpath([Path 'BoundaryValidator/']);
% addpath([Path 'ColorDeconvolution/']);
% addpath([Path 'ColorNormalization/']);
% addpath([Path 'MatOpenSlide/']);

%parse filename
Slashes = strfind(InputFile, '/');
Dots = strfind(InputFile, '.');
if ~isempty(Slashes)
    SlideName = InputFile(Slashes(end)+1:Dots(end)-1);
else
    SlideName = InputFile(1:Dots(end)-1);    
end
OutputFilename = [OutputFolder SlideName '.features.mat'];

%check if file exists
if(~exist(OutputFilename, 'file'))
    
    %check if slide can be opened
    Valid = openslide_can_open(InputFile);
    
    %slide is a valid file
    if(Valid)
        
        %generate schedule for desired magnification
        [Level, Scale, Tout, Factor, X, Y, dX, dY] = ...
            TilingSchedule(InputFile, Objective, T);
        
        %update console
        fprintf('Processing image %s, at magnification %d. Resizing factor %d.\n',...
            InputFile, Objective, Factor);
        
        %initialize containers to capture output from tile analysis
        Features = cell(length(X), 1);
        bX = cell(length(X), 1);
        bY = cell(length(X), 1);
        cX = cell(length(X), 1);
        cY = cell(length(X), 1);
        Elapsed = zeros(length(X), 1);
        
        %process each tilev
        for i = 1:length(X)
            
            %start timer
            tStart = tic;
            
            %read in tile
            I = openslide_read_regions(InputFile, Level, X(i), Y(i),...
                    Tout, Tout);
            
            %update console
            fprintf('\tProcessing tile %06.0f.%06.0f, %d of %d ',...
                    dX(i), dY(i), i, length(X));
            
            %resize if necessary
            if(Factor ~= 1)
                I = imresize(I{1}, Factor, 'bilinear');
            else
                I = I{1};
            end
            
            %normalize color
            RGB = Reinhard(I, TargetMu, TargetSigma, W);
            
            %foreground/background segmentation
            Foreground = NuclearMask(RGB);
            
            %proceed if nuclei are present
            if(sum(Foreground(:)) > 0)
                
                %individual cell segmentation
                [Label, TilebX, TilebY] = NuclearSplitting(Foreground);
                
                %if labeled objects exist
                if(~isempty(TilebX))
                    
                    %feature extraction
                    [TileFeatures, Names, TilecX, TilecY] = ...
                        FeatureExtraction(Label, RGB);
                    
                    %continue processing if objects were located
                    if(~isempty(TileFeatures))
                        
                        %add scan and analysis magnifications to features
                        ScanMag = repmat(Objective/Scale,...
                                [size(TileFeatures,1) 1]);
                        AnalysisMag = repmat(Objective,...
                                [size(TileFeatures, 1) 1]);
                        Features{i} = [ScanMag AnalysisMag TileFeatures];
                        Names = ['ScanMag', 'AnalysisMag', Names];
                        
                        %place boundaries, centroids in global frame,
                        %correct for resizing
                        cX{i} = TilecX + dX(i);
                        cY{i}= TilecY + dY(i);
                        bX{i} = cellfun(@(x) x + dX(i), TilebX,...
                            'UniformOutput', false);
                        bY{i} = cellfun(@(x) x + dY(i), TilebY,...
                            'UniformOutput', false);
                        
                        %stop timer
                        Elapsed(i) = toc(tStart);
                        
                        %write boundary visualization to disk
                        if GenerateImages
                            Mask = bwperim(Label > 0, 4);
                            R = RGB(:,:,1); R(Mask) = 0;
                            G = RGB(:,:,2); G(Mask) = 255;
                            B = RGB(:,:,3); B(Mask) = 0;
                            I = cat(3,R,G,B);
                            Xstr = sprintf('%06.0f', dX(i));
                            Ystr = sprintf('%06.0f', dY(i));
                            imwrite(I, [OutputFolder SlideName '.' ...
                                Xstr '.' Ystr '.jpg']);
                        end
                        
                    end
                end
            end
            
            %update console
            fprintf('%g seconds.\n', toc(tStart));
            
        end
        
        %collapse containers
        Features = cat(1, Features{:});
        bX = cat(1, bX{:});
        bY = cat(1, bY{:});
        cX = cat(1, cX{:});
        cY = cat(1, cY{:});
        
        %discard objects with NaN features
        [Discard, ~] = find(isnan(Features));
        Features(Discard, :) = [];
        bX(Discard) = [];
        bY(Discard) = [];
        cX(Discard) = [];
        cY(Discard) = [];
        
        %calculate sum, sum-squared for normalization
        Sum = nansum(Features, 1);
        SumSq = nansum(Features.^2, 1);
        N = size(Features, 1);
        Mins = min(Features, [], 1);
        Maxes = max(Features, [], 1);
        
        %merge colinear points on boundaries
        for k = 1:length(bX)
            [bX{k}, bY{k}] = MergeColinear(bX{k}, bY{k});
        end
        
        %write results to disk
        save([OutputFolder SlideName '.features.mat'],...
            'bX', 'bY', 'Features', 'Names', 'cX', 'cY', 'Elapsed',...
            'Sum', 'SumSq', 'N', 'Mins', 'Maxes', '-v7.3');
        
        %generate database txt file
        if GenerateText
            SegmentationReport([OutputFolder SlideName '.seg.txt'],...
                SlideName, cX, cY, Features, Names, bX, bY);
        end
        
        %update console
        fprintf('Slide completed - %g seconds.\n', sum(Elapsed));
        
    else
        %display error
        error(['Cannot open slide ' InputFile]);
    end
    
end
