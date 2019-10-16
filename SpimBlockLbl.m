classdef SpimBlockLbl < handle %class of single cell processing of a whole cornea at a single timepoint
    properties
        PartitionNameList
        Frame
        pth
        ImageDims
        tForms
        channelNames
        tile
        
        Centroids
        Intensities
        Int90Prctile
        localEntropyScore
        useInMerge
        NuclearChannel
        
        num
        
        %TopoZ
        %Properties related to tracking
        %Link12Mat
        %Link21Mat
        
    end
    
    properties (Transient = true)
        verbose = true;
    end
    
    properties (Dependent = true)
        filenames
    end
    
    
    
    methods
        
        %         function epiRank = epiRank(W)
        %             [h,x] = hist(W.epiScore,250);
        %             dataCDF = cumsum(h)/sum(h);
        %             epiRank = interp1(x,dataCDF,W.epiScore);
        %             epiRank(isnan(W.epiScore))=NaN;
        %         end
        
        function Blocklbl = SpimBlockLbl(fpath, partitionInfo,timepoint,tile,NuclearChannel)%constructor
            channelNames = partitionInfo.channelNames;
            %Blocklbl = SpimBlockLbl;
            Blocklbl.pth = fpath;
            Blocklbl.Frame = timepoint;
            Blocklbl.tile = tile;
            Blocklbl.NuclearChannel = NuclearChannel;
            %find partitions that correspond to this position and timepoint
            partNameList = partitionInfo.name(find((partitionInfo.timepoint==timepoint).*(partitionInfo.tile==tile)));
            %make sure N files = N channels
            assert(numel(partNameList) == numel(partitionInfo.channelNames));
            Blocklbl.PartitionNameList = partNameList;
            
            
            %Load and process affine transforms
            %each block has a set of corresponding affinetransformations.
            %the transformations should be applied in reverse order (i.e.
            % 5->4->3->2->1). Affine transformations are associative but
            % not commutative, so order must be maintained.
            
            %get all transforms (cell array, 1xnchannels)
            
            tForms = partitionInfo.tForms(find((partitionInfo.timepoint==timepoint).*(partitionInfo.tile==tile)));
            
            
            %get per-channel tranforms
            ChT = cell(numel(tForms),1);
            for j=1:numel(tForms)
                ChT{j} = eye(4);
                for i=1:numel(tForms{j})
                    ChT{j} = ChT{j}*reshape([str2double(strsplit(tForms{j}{i}.affine.Text)), 0 0 0 1]',4,4)';
                end
            end
            
            
            
            
            %Load channels
            options = struct('level', 1,'waitbar',0);
            %Channels = cell(numel(partNameList),1);
            Ints = cell(numel(partNameList),1);
            Ints90P = cell(numel(partNameList),1);
            
            
            %Start calling cells, look at Nuc channel
            indChNuc = find(strcmp(channelNames,NuclearChannel));
            %Data = Channels{indChNuc};
            try
                Data = loadBigDataViewerFormat([fpath filesep partNameList{indChNuc}],options);
                Data = squeeze(single(Data)/2^12);
            catch
                Data = cell(numel(channelNames),1);
                warning('Missing tile on timepoint %d tile %d', timepoint, tile)
            end
            
            if ~isempty(Data)
            localEntropy = 0;%blockproc3(Data,[50 50, 50],@(x) simpEntropy(x));
            Blocklbl.ImageDims = size(Data);
            
            % Gaussian smoothen and find regional max
            
            imgG = fastGauss3D(Data,[4,4,1]);
            %imghMaxima = imhmax(imgG, 0.00001);
            %imgGRegMax = imregionalmax(imgG);
            
            
            localMax = dlmread(fullfile(fpath, [partitionInfo.beadsFile{find((partitionInfo.timepoint==timepoint).*(partitionInfo.tile==tile).*(partitionInfo.channel==indChNuc))}]),'\t',2,1);
            localMax((localMax(:,2)<1),1)=1;
            localMax((localMax(:,1)<1),2)=1;
            localMax((localMax(:,3)<1),3)=1;
            localMax((localMax(:,2)>size(Data,1)),1)=size(Data,1);
            localMax((localMax(:,1)>size(Data,2)),2)=size(Data,2);
            localMax((localMax(:,3)>size(Data,3)),3)=size(Data,3);
            
            linearInd = sub2ind(size(Data), round(localMax(:,2)), round(localMax(:,1)), round(localMax(:,3)));
            imgGRegMax = zeros(size(Data));
            imgGRegMax(linearInd)=1;
            % Expand regional max to a sphere of radius=4
            [xx,yy,zz] = ndgrid(-3:3);
            nhood = sqrt(xx.^2 + yy.^2+zz.^2) <= 3.0;
            ExpandedPeaks = imdilate(imgGRegMax,nhood);
            clear imgGRegMax  nhood
            % Filter peaks with low signal using triangle thresholding
            somePix = imgG(logical(ExpandedPeaks));
            %[h,x]=hist(log(somePix),100);
            %ThreshVal = exp(triangleThresh(h,x));
            ThreshVal = exp(meanThresh(log(somePix)));
            %ThreshVal = (meanThresh((somePix)));
            nuclearSeeds = ExpandedPeaks.*(imgG>ThreshVal);
            %nuclearSeedsSignal = Data.*nuclearSeeds;
            clear ExpandedPeaks imgG somePix
            % Find connected components
            CC = bwconncomp(nuclearSeeds);
            clear nuclearSeeds
            % get region props in each channel
            S = regionprops(CC,Data,'MeanIntensity','WeightedCentroid','Area','PixelValues');
            
            clear Data
            
            % filter CC keep only cells with volume in some range
            Areas = cat(1, S.Area);
            Intensities = cat(1, S.MeanIntensity).*Areas;
            Int90Prctile = arrayfun(@(x) prctile(x.PixelValues(x.PixelValues~=0),90),S);
            Centroids = cat(1,S.WeightedCentroid);
            %remove specks
            J = Areas>=120;
            Areas = Areas(J);
            Intensities = Intensities(J);
            Int90Prctile = Int90Prctile(J);
            Centroids = Centroids(J,:);
            
            %Add to Int cell array
            Ints{indChNuc}=Intensities;
            Ints90P{indChNuc}= Int90Prctile;
            
            
            
            %S = regionprops(CC,localEntropy,'MeanIntensity');
            %localEntropyScore = cat(1, S.MeanIntensity);
            %localEntropyScore = localEntropyScore(J);
            localEntropyScore = 0;
            
            
            %Now the other channels
            indOtherChannels = find(~strcmp(channelNames,NuclearChannel));
            for i=indOtherChannels
                
                try
                    Data = loadBigDataViewerFormat([fpath filesep partNameList{i}],options);
                catch
                    Ints{i}=nan(nnz(J),1);
                    Ints90P{i}= nan(nnz(J),1);
                    continue
                end
                Data = squeeze(single(Data)/2^12);
                %RI = imref3d(size(Data));
                
                
                
                %%%
                %we're gonna try a different strategy. In stead of
                %transforming the data, we'll move the points we found
                %before.
                %[i,j,k]=ind2sub(CC.ImageSize,CC.PixelIdxList{1})
                %%%%%
                %This is all wrong
                %                 Data = imwarp(Data,RI,affine3d(tform'));
                %
                %                 %pad missing bits with zeros
                %                 wXYZ = size(Data);
                %                 dwXYZ = CC.ImageSize-wXYZ;
                %                 padSize = dwXYZ.*(dwXYZ>0);
                %                 Data = padarray(Data,padSize,0,'post');
                %                 %crop extra bits
                %                 Data = Data(1:CC.ImageSize(1),1:CC.ImageSize(2),1:CC.ImageSize(3));
                %                 %%%%%
                tform = inv(ChT{indChNuc})*ChT{i};
                Tform3D = affine3d(tform');
                CC0 = CC;
                
                for ind=1:CC.NumObjects
                    [ix,iy,iz]=ind2sub(CC0.ImageSize,CC.PixelIdxList{ind});
                    [ii,jj,kk]=Tform3D.transformPointsInverse(ix,iy,iz);
                    ii=round(ii);
                    jj=round(jj);
                    kk=round(kk);
                    Jtokeep = logical((ii>0).*(jj>0).*(kk>0).*(ii<=CC0.ImageSize(1)).*(jj<=CC0.ImageSize(2)).*(kk<=CC0.ImageSize(3)));
                    inds = unique(sub2ind(CC0.ImageSize,ii(Jtokeep),jj(Jtokeep),kk(Jtokeep)));
                    CC0.PixelIdxList{ind} = inds;
                end
                
                
                S = regionprops(CC0,Data,'MeanIntensity','PixelValues');
                Intensities = cat(1, S.MeanIntensity);
                Int90Prctile = arrayfun(@(x) prctile(x.PixelValues(x.PixelValues~=0),90),S);
                %J is still the index of no specky cells
                Ints{i}=Intensities(J);
                Ints90P{i}= Int90Prctile(J);
                clear Data;
            end
            clear CC CC0;
            
            
            
            %Apply global transformations to centroids!
            tempCents = (ChT{indChNuc}*[Centroids ones(1,size(Centroids,1))']')';
            Centroids = tempCents(:,1:3);
            else
                Centroids = [];
                Ints = [];
                Ints90P = [];
                localEntropyScore = [];
            end
            
            Blocklbl.Centroids = Centroids;
            Blocklbl.Intensities = Ints;
            Blocklbl.Int90Prctile = Ints90P;
            Blocklbl.localEntropyScore = localEntropyScore;
            Blocklbl.useInMerge = true(size(Centroids,1),1);
            
            Blocklbl.tForms = tForms;
            Blocklbl.channelNames = channelNames;
            Blocklbl.num = size(Centroids,1);
            %Blocklbl.DistFromManifold = DistFromManifold;
            %Blocklbl.TopoZ = TopoZ;
            %Blocklbl.GridZ = zz;
            %Blocklbl.CC = CC;
            tile
        end
        
        function bb = boundingBox(W)
            if W.num>1
                bb = [min(W.Centroids)' max(W.Centroids)'];
            else
                bb=[0 0 0; 0 0 0]';
            end
        end
        
        
        function E = simpEntropy(J1)
            h = histcounts(J1,'Normalization','probability');
            
            E = -nansum(h.*log(h));
        end
        
        
        function mem = freeMem(W)
            [~,out]=system('vmstat -s -S M | grep "free memory"');
            
            mem=sscanf(out,'%f  free memory');
        end
        
        
        function h = scatter3(W,varargin)
            
            channelToShow = ParseInputs('channel', [], varargin);
            channelDenom = ParseInputs('channel2', [], varargin);
            
            ratio = ParseInputs('ratio', false, varargin);
            Jplot = ParseInputs('range', 1:W.num, varargin);
            dz = ParseInputs('dz', 0, varargin);
            
            %first, decide which points to plot
            if strcmp(Jplot,'useInMerge')
                Jplot = find(W.useInMerge);
            end
            if any(Jplot>W.num)
                Jplot = Jplot(Jplot<=W.num);
                disp('Range out of bounds, drawing all valid points.')
            end
            
            %next, decide on channel
            
            if isempty(channelToShow);
                tzeva = viridis(numel(Jplot));
            elseif ~ratio
                indChNuc = find(strcmp(W.channelNames,channelToShow));
                %tzeva = viridis(length(W.Centroids));
                tzeva = W.Intensities{indChNuc};
                tzeva = tzeva(Jplot);
            else
                if isempty(channelDenom)
                    error('need 2 channels for ratio')
                end
                indChNuc = find(strcmp(W.channelNames,channelToShow));
                indChDenom = find(strcmp(W.channelNames,channelDenom));
                
                %tzeva = viridis(length(W.Centroids));
                tzeva = W.Intensities{indChNuc}./W.Intensities{indChDenom};
                tzeva = tzeva(Jplot);
            end
            
            h = scatter3(W.Centroids(Jplot,1),W.Centroids(Jplot,2),dz-W.Centroids(Jplot,3),[],tzeva);
            
        end
        
        
        
        
        function scatter(W,varargin)
            tzeva = viridis(length(W.Centroids));
            if nargin==1
                scatter(W.Centroids(:,1),W.Centroids(:,2),[],tzeva);
            else
                range = varargin{1};
                if strcmp(range,'epi')
                    tzeva = viridis(length(W.Jepi));
                    if nargin>=3;
                        range = varargin{2};
                        if any(~ismember(range,W.Jepi))
                            range = range(ismember(range,W.Jepi));
                            if ~isempty(range)
                                disp('some of the range points are not in the epithelium, drawing the ones that are.')
                            else
                                disp('all of the range points are not in the epithelium, drawing entire epithelium.')
                                range = W.Jepi;
                            end
                        end
                    else
                        range = W.Jepi;
                    end
                    if any(range>max(W.Jepi))
                        range = W.Jepi;
                        disp('range out of bounds, drawing all points.')
                        scatter(W.Centroids(range,1),W.Centroids(range,2),[],tzeva((ismember(W.Jepi,range)),:));
                    else
                        scatter(W.Centroids(range,1),W.Centroids(range,2),[],tzeva((ismember(W.Jepi,range)),:));
                    end;
                else
                    if any(range>W.num)
                        range = range(1):W.num;
                        disp('range out of bounds, drawing all points up to the total # of cells.')
                        scatter(W.Centroids(range,1),W.Centroids(range,2),[],tzeva(range,:));
                    else
                        scatter(W.Centroids(range,1),W.Centroids(range,2),[],tzeva(range,:));
                    end;
                end
            end
            %scatter(Centroids(:,1),Centroids(:,2),[],parula(length(Centroids)));
            %set(gca,'xlim',[-200 3000],'ylim',[-200 2500])
            %hold on
            shg
            
        end
        
        
        
        
        
        function stkshow(W)
            MD=Metadata(W.pth);
            
            Data =  stkread(MD,'Channel','DeepBlue', 'flatfieldcorrection', false, 'frame', W.Frame, 'Position', W.PosName,'register',false);
            RChannel=Data;
            GChannel=Data;
            BChannel=Data;
            [h,x] = hist(log(datasample(Data(:),min(1000000,numel(Data(:))))),1000);
            maxC = exp(x(find(cumsum(h)./sum(h)>0.995,1,'first')));
            cmap = parula(W.CC.NumObjects)*maxC*1.5;%scale colormap so that data and centroids are all visible
            
            for i=1:W.CC.NumObjects
                indexToChange = W.CC.PixelIdxList{i};
                RChannel(indexToChange)=cmap(i,1); %Replace actual pixel values with relevant RGB value
                GChannel(indexToChange)=cmap(i,2);
                BChannel(indexToChange)=cmap(i,3);
                i
            end
            RGB = cat(3,RChannel,GChannel,BChannel);%combine into a single RGB image and show. stkshow is sloooooooooooowing me down.
            stkshow(RGB);
            MIJ.selectWindow('RGB');
            MIJ.run('Stack to Hyperstack...', ['order=xyzct channels=3 slices=' num2str(size(Data,3)) ' frames=1 display=Composite']);
            %stkshow(RChannel)
            %stkshow(GChannel)
            %stkshow(BChannel)
            %MIJ.run('Merge Channels...', 'c1=RChannel c2=GChannel c3=BChannel create');
        end
        
        
        function scattershow(W,varargin)
            MD=Metadata(W.pth);
            
            Data =  stkread(MD,'Channel','DeepBlue', 'flatfieldcorrection', false, 'frame', W.Frame, 'Position', W.PosName,'register',false);
            RChannel=zeros(size(Data));
            GChannel=zeros(size(Data));
            BChannel=zeros(size(Data));
            [h,x] = hist(log(datasample(Data(:),min(1000000,numel(Data(:))))),1000);
            maxC = exp(x(find(cumsum(h)./sum(h)>0.995,1,'first')));
            if strcmp(varargin{1},'epi')
                range = W.Jepi'
            else
                range = 1:W.CC.NumObjects;
            end
            cmap = viridis(numel(range))*maxC*1.5;%scale colormap so that data and centroids are all visible
            
            for i=1:numel(range)
                indexToChange = W.CC.PixelIdxList{range(i)};
                RChannel(indexToChange)=cmap(i,1); %Replace actual pixel values with relevant RGB value
                GChannel(indexToChange)=cmap(i,2);
                BChannel(indexToChange)=cmap(i,3);
                i
            end
            RGB = cat(3,RChannel,GChannel,BChannel);%combine into a single RGB image and show. stkshow is sloooooooooooowing me down.
            stkshow(RGB);
            MIJ.selectWindow('RGB');
            MIJ.run('Stack to Hyperstack...', ['order=xyzct channels=3 slices=' num2str(size(Data,3)) ' frames=1 display=Composite']);
        end
        
    end
end
