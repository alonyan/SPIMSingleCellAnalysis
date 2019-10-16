classdef MultiColorSPIMTimepoint
    properties
        pth
        BlockLbls
        timepoint
        NuclearChannel
        channelNames
        Centroids
        Intensities
        Int90Prctile
        num
        Link12
        Link21
    end
    
    properties (Transient = true)
        verbose = true;
    end
    
    properties (Dependent = true)
        filenames
    end
    
    properties (Hidden = true)
        notBeads
    end
    
    
    methods
        %constructor
        function W = MultiColorSPIMTimepoint(fpath,partitionInfo, timepoint,NuclearChannel,varargin)
            override = ParseInputs('override', false, varargin);
            foldname = fullfile(fpath,'SingleCellSPIMTimepoints/');
            filename = sprintf('MultiColorSPIMTimepoint_%03d.mat',timepoint);
            filename = fullfile(foldname,filename);
            if exist(filename,'file') && ~override
                s=load(filename);
                W = s.W;                              
                %% we're done.
                return
            else              
                tilesToLoad = unique(partitionInfo.tile((partitionInfo.timepoint)==timepoint));
                BlockLbls = cell(numel(tilesToLoad),1);
                parfor i=1:numel(tilesToLoad)
                    %this sets the # of cores per worker to 4!
                    warning('off')
                    maxNumCompThreads(4);
                    workersCompThreads(i) = maxNumCompThreads;
                    warning('on')

                    BlockLbls{i} = SpimBlockLbl(fpath, partitionInfo,timepoint,tilesToLoad(i),NuclearChannel);
                end
                
                %             spmd
                %                 origLimit = numel(tilesToLoad);
                %                 loopLimit = numlabs * ceil(origLimit/numlabs);
                %                 for i = labindex:numlabs:loopLimit
                %                     %labBarrier; % synchronise all workers
                %                     if i <= origLimit
                %                      %   pause(30*(labindex-1)); % timing offset
                %                         BlockLbls{i+1} = SpimBlockLbl(fpath, partitionInfo,timepoint,tilesToLoad(i),NuclearChannel);
                %                     end
                %                 end
                %             end
                
                
                W.BlockLbls = BlockLbls;
                W.channelNames = W.BlockLbls{1}.channelNames;
                W.pth = fpath;
                W.timepoint = timepoint;
                W.NuclearChannel = NuclearChannel;
                
                %combine all the centroids from all the tiles
                W.mergeTiles;
                W.save;
            end
        end
        
        
        
        
        function save(W)
            foldname = fullfile(W.pth,'SingleCellSPIMTimepoints');
            if ~exist(foldname,'dir')
                mkdir(foldname)
            end
            filename = sprintf('MultiColorSPIMTimepoint_%03d',W.timepoint);
            filename = fullfile(foldname,filename);
            save(filename, 'W');
        end
        
        function num = get.num(W)
            num = size(W.Centroids,1);
        end
        
        
        function Cents = get.Centroids(W)
            Cents = cellfun(@(x) x.Centroids(logical(x.useInMerge),:), W.BlockLbls,'uniformOutput',false);
            Cents = cat(1,Cents{:});
            notBeads = W.findNotBeads(Cents);
            Cents = Cents(notBeads,:);
        end
        
        function IntsOut = get.Intensities(W)
            Ints = cellfun(@(y) cellfun(@(x) x(logical(y.useInMerge)), y.Intensities, 'uniformoutput', false), W.BlockLbls, 'uniformoutput', false);
            Ints = cat(1,Ints{:});
            Ints = reshape(Ints,numel(W.channelNames),[]);
            IntsOut = cell(numel(W.channelNames),1);
            for i=1:numel(W.channelNames)
                IntsOut{i} = cat(1,Ints{i,:});
                IntsOut{i} = IntsOut{i}(W.notBeads);
            end
            clear Ints;
        end
        
        function IntsOut = get.Int90Prctile(W)
            Ints = cellfun(@(y) cellfun(@(x) x(logical(y.useInMerge)), y.Int90Prctile, 'uniformoutput', false), W.BlockLbls, 'uniformoutput', false);
            Ints = cat(1,Ints{:});
            Ints = reshape(Ints,numel(W.channelNames),[]);
            IntsOut = cell(numel(W.channelNames),1);
            for i=1:numel(W.channelNames)
                IntsOut{i} = cat(1,Ints{i,:});
                IntsOut{i} = IntsOut{i}(W.notBeads);
            end
            clear Ints;
        end
        
        
        
        
        
        
        
        
        
        
        function mergeTiles(W)
            nTiles = numel(W.BlockLbls);
            for i=1 : nTiles
                i
                for j= i+1 : nTiles
                    if W.blocksoverlap(i,j)
                        W.matchBlocks(i,j);
                    end
                end
                
            end
        end
        
        
        %function that checks whether 2 blocks overlap so we can avoid
        %unnecessary calculations
        function overlapqustionmark = blocksoverlap(W, i, j)
            bounds1 = W.BlockLbls{i}.boundingBox;
            bounds2 = W.BlockLbls{j}.boundingBox;
            dim=2;
            
            min1 = min(bounds1,[],dim);
            max1 = max(bounds1,[],dim);
            min2 = min(bounds2,[],dim);
            max2 = max(bounds2,[],dim);
            d = max(0, min([max1 , max2],[],dim)-max([min1 , min2],[],dim));
            overlapqustionmark = prod(d)>0;
        end
        
        %function that looks at 2 blocks, finds matched points, and labels
        %which one to use in the merge
        function matchBlocks(W,i,j)
            %get 2 blocks
            block1 = W.BlockLbls{i};
            block2 = W.BlockLbls{j};
            
            %index of nuclear channel
            indChNuc = find(strcmp(block1.channelNames,block1.NuclearChannel));
            
            
            searchRadius=50;
            
            %find 5 nearest neighbors in block2 of points in block1
            [n,d]=knnsearch(block2.Centroids,block1.Centroids,'k',5);
            
            %find the mean distance
            d=median(d,2);
            
            %find indexes of centroids in block1 that have close neighbors in block 2.
            indD = find(d<=searchRadius);
            %find the indexes of those neighbors and calculate their mean intensity
            a = n(indD,:);
            a = mat2cell(a,ones(size(a,1),1),size(a,2));
            medianIntensityB = cellfun(@(x) mean(block2.Intensities{indChNuc}(x)), a);
            
            %intensity of matched points in block1
            IntensityA = block1.Intensities{indChNuc}(indD);
            
            block1.useInMerge(indD(IntensityA<=medianIntensityB)) = false;
            
            %repeat for other direction
            %find 5 nearest neighbors in block2 of points in block1
            [n,d]=knnsearch(block1.Centroids,block2.Centroids,'k',5);
            
            %find the mean distance
            d=median(d,2);
            
            %find indexes of centroids in block1 that have close neighbors in block 2.
            indD = find(d<=searchRadius);
            %find the indexes of those neighbors and calculate their mean intensity
            a = n(indD,:);
            a = mat2cell(a,ones(size(a,1),1),size(a,2));
            medianIntensityB = cellfun(@(x) mean(block1.Intensities{indChNuc}(x)), a);
            
            %intensity of matched points in block1
            IntensityA = block2.Intensities{indChNuc}(indD);
            
            block2.useInMerge(indD(IntensityA<=medianIntensityB)) = false;
            
            
            
            %Next look for points that exactly match
            %look within a 15 pixel radius for matching points
            searchRadius=15;
            Dists = createDistanceMatrix(block1.Centroids, block2.Centroids);
            costMat = Dists;
            
            costMat(Dists>searchRadius) = 0;
            costMat = sparse(double(costMat));
            [Links12, ~] = lap(costMat,[],[],1);
            
            
            Link12Mat =repmat(Links12(1:length(block1.Centroids(:,1))),1,length(block2.Centroids(:,1)));
            LinkMat1 = meshgrid(1:length(block2.Centroids(:,1)), 1:length(block1.Centroids(:,1)));
            Link12Mat = LinkMat1==Link12Mat;
            clear LinkMat1 costMat;
            [ind1, ind2]=find(Link12Mat);
            
            indChNuc = find(strcmp(block1.channelNames,block1.NuclearChannel));

            %select the points with higher local entropy to be used
            block1.useInMerge(ind1(block1.Intensities{indChNuc}(ind1) < block2.Intensities{indChNuc}(ind2)))=false;
            block2.useInMerge(ind1(block1.Intensities{indChNuc}(ind1) >= block2.Intensities{indChNuc}(ind2)))=false;
            
            W.BlockLbls{i} = block1;
            W.BlockLbls{j} = block2;
        end
        
        function notBeads = get.notBeads(W)
            Cents = cellfun(@(x) x.Centroids(logical(x.useInMerge),:), W.BlockLbls,'uniformOutput',false);
            Cents = cat(1,Cents{:});
            notBeads = W.findNotBeads(Cents);
        end
        function notBeads = findNotBeads(W, Centroids)
            [~,d]=knnsearch(Centroids,Centroids,'k',6);
            notBeads = mean(d(:,2:end),2)<75;       
        end
        
        function D = Density(W,J)
           [~,D] = knnsearch(W.Centroids,W.Centroids,'K', 31); 
           D = 1./mean(D(:,2:end),2);
           if nargin==1
               J=1:numel(D);
           end
           D = D(J);
        end
        
        function h= scatter3(W,varargin)
            
            channelToShow = ParseInputs('channel', [], varargin);
            channelDenom = ParseInputs('channel2', [], varargin);
          
            thresh = ParseInputs('thresh', 0, varargin);            
            ratio = ParseInputs('ratio', false, varargin);
            Jplot = ParseInputs('range', 1:W.num, varargin);
            dz = ParseInputs('dz', 0, varargin);
            
            %which points
            if any(Jplot>W.num)
                Jplot = Jplot(Jplot<=W.num);
                disp('Range out of bounds, drawing all valid points.')
            end
            
            
            %next, decide on channel
            
            if isempty(channelToShow); 
                %tzeva = viridis(numel(Jplot));
                if numel(Jplot)>1
                    tzeva = [viridis(ceil(numel(Jplot)/2)); flipud(viridis(floor(numel(Jplot)/2)))];
                else
                    tzeva=[0 0 1];
                end

            elseif ~ratio
                indChNuc = find(strcmp(W.channelNames,channelToShow));
                %tzeva = viridis(length(W.Centroids));
                tzeva = log(W.Int90Prctile{indChNuc});
                if thresh
                    tzeva = W.Int90Prctile{indChNuc}>thresh;
                end
                tzeva = tzeva(Jplot);

            else
                if isempty(channelDenom)
                    error('need 2 channels for ratio')   
                end
                indChNuc = find(strcmp(W.channelNames,channelToShow));
                indChDenom = find(strcmp(W.channelNames,channelDenom));

                %tzeva = viridis(length(W.Centroids));
                tzeva = W.Int90Prctile{indChNuc}./W.Int90Prctile{indChDenom};
                tzeva = tzeva(Jplot);
            end
            
            
            
            h = scatter3(W.Centroids(Jplot,1),W.Centroids(Jplot,2),dz-W.Centroids(Jplot,3),3,tzeva,'filled');

           
            %Ch = get(gca,'Children');
            %nPoints = numel(Ch.XData);
            %scatter(Centroids(:,1),Centroids(:,2),[],parula(length(Centroids)));
            set(gca,'CameraPositionMode','manual','CameraPosition', [3.9379e+04 -5.3724e+03 3.9412e+03]...
                ,'color','k', 'xcolor','none', 'ycolor','none', 'zcolor','none','xgrid', 'off','ygrid', 'off','zgrid', 'off')
            set(gcf,'color','k')
            set(gca,'xlim',[[4000 8000]],'ylim',[1000 7000],'zlim',[-4000 500])
            hold on
            %shg
            
        end
        
    end
    
end