classdef MulticolorSPIMResults < Results
    properties
        Frames = {};
        SPIMTimepoints
        partitionInfo
        Tracks
        NuclearChannel
    end
    
    properties (Dependent = true)
        TimeVecs
        loadTPsflag
    end
    
    methods (Static)
        function R = loadobj(S)
            R = MulticolorSPIMResults;
            R = reload(R,S);
        end
    end
    
    
    methods
        
        function S = saveobj(R)
            S = toStruct(R);
        end
        
        function R = reload(R,S)
            R = reload@Results(R,S);
            
            if isfield(S,'Frames') % is the fields exists load them, if not,
                % they will be filled with default values
                % effectivly upcasting an object.
                R.Frames = S.Frames;
                %R.Frames = [0:72]';
            end
            if isfield(S,'SPIMTimepoints') % is the fields exists load them, if not,
                % they will be filled with default values
                % effectivly upcasting an object.
                %R.WoundLbl = S.WoundLbl;
                R.SPIMTimepoints = {};
                if R.loadTPsflag
                    SPIMTPs = cell(1,numel(R.Frames));
                    
                    for i=R.Frames'
                        disp(['loading timepoint ' num2str(i)])
                        SPIMTPs{i+1} = R.getSPIMTimepoints(i);
                    end
                    R.SPIMTimepoints = SPIMTPs;
                end
            end
            if isfield(S,'Tracks') % is the fields exists load them, if not,
                % they will be filled with default values
                % effectivly upcasting an object.
                R.Tracks = S.Tracks;
            end
            if isfield(S,'partitionInfo') % is the fields exists load them, if not,
                % they will be filled with default values
                % effectivly upcasting an object.
                R.partitionInfo = S.partitionInfo;
            end
        end
        
        function S = toStruct(R)
            % call the superclass method to start the transition to a
            % struct
            S = toStruct@Results(R);
            
            %S.WoundLbl = R.WoundLbl;
            S.SPIMTimepoints = R.SPIMTimepoints;
            S.Tracks = R.Tracks;
            S.Frames = R.Frames;
            S.partitionInfo = R.partitionInfo;
            %             S.TimeVecs = R.TimeVecs;
        end
        
        function R = MulticolorSPIMResults(pth,reset) %constructor
            if nargin==0
                pth='';
                reset=false;
            end
            if nargin==1
                reset=false;
            end
            R@Results(pth,reset);
        end
        
        %function TimeVecs = get.TimeVecs(R)
        %    TimeVecs = getTimeVecs(R);
        %end
        
        function TimeVecs = get.TimeVecs(R)
            TimeVecs{ix} = unique(R.partitionInfo.timepoint);
        end
        
        function loadTPsflag = get.loadTPsflag(R)
            foldname = fullfile(R.pth,'SingleCellSPIMTimepoints/');
            loadTPsflag = exist(foldname,'dir');
        end
        
        function set.SPIMTimepoints(R, SPIMTPs)
            R.SPIMTimepoints = SPIMTPs;
        end
        
        function SPIMTP = calculateSPIMTimepoints(R,timepoint,varargin)
            SPIMTP = MultiColorSPIMTimepoint(R.pth,R.partitionInfo, timepoint,R.NuclearChannel,varargin{:});
        end
        
        function SPIMTP = getSPIMTimepoints(R,timepoint)
            foldname = fullfile(R.pth,'SingleCellSPIMTimepoints/');
            filename = sprintf('MultiColorSPIMTimepoint_%03d.mat',timepoint);
            filename = fullfile(foldname,filename);
            if exist(filename,'file')
                SPIMTP = MultiColorSPIMTimepoint(R.pth,R.partitionInfo, timepoint,R.NuclearChannel);
            else
                SPIMTP = {};
                warning('No MultiColorSPIMTimepoint found on timepoint %d' , timepoint)
            end
        end
        
        function Link(R,trckChnl) %need to revise for bigger cell#
            %% Now, we'll see how well we can track between adjacent frames (lap)
            
            SPIMTP = R.SPIMTimepoints;
            
            %init assignment matrices
            Link12MatCell = {};
            Link21MatCell = {};
            
            searchRadius = 100;
            maxAmpRatio = 6;
            epsilon = 10^-4;%prevent cost==0
            
            inds = find(cellfun(@(x) ~isempty(x), SPIMTP));
            TrackChInd = find(strcmp(trckChnl,SPIMTP{inds(1)}.channelNames));
            
            for i=inds(1:end-1);
                i
                SPIMTPi = SPIMTP{i};
                SPIMTPip1 = SPIMTP{i+1};
                
                
                ND = SPIMTPi.num;
                NA = SPIMTPip1.num;
                CentsD = SPIMTPi.Centroids;
                CentsA = SPIMTPip1.Centroids;
                IntD = SPIMTPi.Intensities{TrackChInd};
                IntA = SPIMTPip1.Intensities{TrackChInd};
                
                
                [ids , dist] = rangesearch(CentsA,CentsD,searchRadius);
                
                indx1 = arrayfun(@(x) repmat(x,1,numel(ids{x})), 1:numel(ids),'UniformOutput', false);
                indx1 = cat(2,indx1{:});
                indx2 = cat(2,ids{:});
                Dists = cat(2,dist{:});
                ampRatio = max(IntD(indx1),IntA(indx2))./min(IntD(indx1),IntA(indx2));
                
                costs = double(Dists.*(log2(ampRatio')+epsilon));
                
                
                
                
                costMat = sparse(indx1,indx2,costs,ND,NA);
                
                costMat(find(isnan(costMat)))=0;
                
                clearvars costs indx1 indx2 SPIMTPi SPIMTPip1 CentsD CentsA IntD IntA
                
                [Links12, Links21] = lap(costMat,[],[],1);
                
                Link12 = int32(Links12(1:ND)).*int32(Links12(1:ND)<=NA);
                Link21 = int32(Links21(1:NA)).*int32(Links21(1:NA)<=ND);
                
                SPIMTP{i}.Link12  = Link12;
                SPIMTP{i}.Link21 = Link21;
                SPIMTP{i}.save;
            end
            
            
            R.SPIMTimepoints = SPIMTP;
        end
        
        function closeGaps(R,trckChnl)
            SPIMTP = R.SPIMTimepoints;
            indtrckChnl = find(strcmp(trckChnl,SPIMTP{1}.channelNames));
            
            inds = find(cellfun(@(x) ~isempty(x), SPIMTP));
            %% Book keeping: Make all tracks fragmants
            numFeatures = cellfun(@(x) x.num, SPIMTP(inds))';
            trackedFeatureIndx = (1:numFeatures(1))';
            numFrames = numel(inds);
            
            %initialize auxiliary matrices for storing information related to tracks
            %fragments
            numTracksWorstCase = round(sum(numFeatures)/10); %arbitrary large number
            trackedFeatureIndxAux = zeros(numTracksWorstCase,numFrames);
            rowEnd = numTracksWorstCase; %We'll fill this from the bottom up
            
            for i=1:numFrames-1
                %get indices of features in 2nd frame that are connected to features in 1st frame
                %indx1C - indexes in frame 1 that are linked to frame 2
                %indx2C - indexes in frame 2 that are linked to indx1C in frame 1
                numFeaturesFrame1 = numFeatures(i);
                numFeaturesFrame2 = numFeatures(i+1);
                %[indx2C,indx1C] = find(SPIMTP{i}.Link21Mat);
                indx2C = find(SPIMTP{i}.Link21);
                indx1C = SPIMTP{i}.Link21(indx2C);
                %%
                %find existing tracks that are not connected to features in 2nd frame
                numExistTracks = size(trackedFeatureIndx,1);
                indx1U = setdiff(1:numExistTracks,indx1C); %features in 1 not connected to 2
                numRows = length(indx1U);
                %%
                %determine where to store these tracks in auxiliary matrix
                %extend auxiliary matrices if necessary
                rowStart = rowEnd - numRows + 1;
                if rowStart <= 1
                    trackedFeatureIndxAux = [zeros(numTracksWorstCase,numFrames); ...
                        trackedFeatureIndxAux];
                    rowEnd = rowEnd + numTracksWorstCase;
                    rowStart = rowStart + numTracksWorstCase;
                end
                
                %% move rows of tracks that are not connected to points in
                %2nd frame to auxilary matrix
                trackedFeatureIndxAux(rowStart:rowEnd,1:i) = trackedFeatureIndx(indx1U,:);
                
                %%
                %assign space for new connectivity matrix
                tmp = zeros(numFeaturesFrame2,i+1);
                %fill in the feature numbers in 2nd frame
                tmp(1:numFeaturesFrame2,i+1) = (1:numFeaturesFrame2)';
                %shuffle existing tracks to get the correct connectivity with 2nd frame
                tmp(indx2C,1:i) = trackedFeatureIndx(indx1C,:);
                %update the connectivity matrix "trackedFeatureIndx"
                trackedFeatureIndx = tmp;
                
                %update rowEnd to indicate until which row the auxiliary
                %matrices are ampty
                rowEnd = rowStart - 1;
            end
            %
            %add information from last frame to auxiliary matrices
            numRows = size(trackedFeatureIndx,1);
            rowStart = rowEnd - numRows + 1;
            if rowStart <= 1
                trackedFeatureIndxAux = [zeros(numRows,numFrames); ...
                    trackedFeatureIndxAux];
                rowEnd = rowEnd + numRows;
                rowStart = rowStart + numRows;
            end
            trackedFeatureIndxAux(rowStart:rowEnd,:) = trackedFeatureIndx;
            
            %remove all empty rows
            trackedFeatureIndx = trackedFeatureIndxAux(rowStart:end,:);
            clear trackedFeatureIndxAux
            
            % get total number of tracks
            numTracks = size(trackedFeatureIndx,1);
            
            %find the frame where each track begins and then sort the vector
            frameStart = zeros(numTracks,1);
            for i=1:numTracks
                frameStart(i) = find((trackedFeatureIndx(i,:)~=0),1,'first');
            end
            [frameStart,indx] = sort(frameStart);
            
            %rearrange "trackedFeatureIndx" such that tracks are sorted in ascending order by their
            %starting point. Note that this ends up also arranging tracks starting at the
            %same time in descending order from longest to shortest.
            trackedFeatureIndx = trackedFeatureIndx(indx,:);
            
            
            %% Filter short tracks
            movieInfo = [];
            for i=inds
                n = SPIMTP{i}.num;
                movieInfo(i).xCoord = [SPIMTP{i}.Centroids(:,1) zeros(n,1)];
                movieInfo(i).yCoord = [SPIMTP{i}.Centroids(:,2) zeros(n,1)];
                movieInfo(i).zCoord = [SPIMTP{i}.Centroids(:,3) zeros(n,1)];
                movieInfo(i).amp = [SPIMTP{i}.Int90Prctile{indtrckChnl}  zeros(n,1)];
                movieInfo(i).num = SPIMTP{i}.num;
            end
            probDim = 3;
            trackedFeatureInfo = coordAmpMatFromIndicesSparse(trackedFeatureIndx,movieInfo,...
                numFrames,probDim);
            trackSEL = getTrackSEL(trackedFeatureInfo);
            
            minTrackLen=3;
            
            %remove tracks whose length is less than minTrackLen
            indxKeep = find(trackSEL(:,3) >= minTrackLen);
            trackSEL = trackSEL(indxKeep,:);
            trackedFeatureIndx = trackedFeatureIndx(indxKeep,:);
            trackedFeatureInfo = trackedFeatureInfo(indxKeep,:);
            numTracks = size(trackSEL,1)
            clear movieInfo;
            
            
            
            
            %% Find all possible links based on thresholds
            maxTimeJump = 3;
            maxStep = 75;
            trackStartTime = trackSEL(:,1);
            trackEndTime   = trackSEL(:,2);
            %CentroidsStarts = [];
            %CentroidsEnds = [];
            
            %             for ind=1:numTracks
            %                 CentroidsStarts = [CentroidsStarts; full(trackedFeatureInfo(ind,8*(trackStartTime(ind)-1)+1:8*(trackStartTime(ind)-1)+3))];
            %                 CentroidsEnds = [CentroidsEnds; full(trackedFeatureInfo(ind,8*(trackStartTime(ind)-1)+1:8*(trackStartTime(ind)-1)+3))];
            %             end
            CentroidsStarts = cell2mat(arrayfun(@(a) full(trackedFeatureInfo(a,8*(trackStartTime(a)-1)+1:8*(trackStartTime(a)-1)+3))',1:numTracks,'UniformOutput',false))';
            CentroidsEnds = cell2mat(arrayfun(@(a) full(trackedFeatureInfo(a,8*(trackEndTime(a)-1)+1:8*(trackEndTime(a)-1)+3))',1:numTracks,'UniformOutput',false))';
            
            AmpStarts = cell2mat(arrayfun(@(a) full(trackedFeatureInfo(a,8*(trackStartTime(a)-1)+4))',1:numTracks,'UniformOutput',false))';
            AmpEnds = cell2mat(arrayfun(@(a) full(trackedFeatureInfo(a,8*(trackEndTime(a)-1)+4))',1:numTracks,'UniformOutput',false))';
            
            
            searchRadius = 100;
            
            [ids , dist] = rangesearch(CentroidsEnds, CentroidsStarts,searchRadius);
            %
            indx1 = arrayfun(@(x) repmat(x,1,numel(ids{x})), 1:numel(ids),'UniformOutput', false);
            indx1 = cat(2,indx1{:});
            indx2 = cat(2,ids{:});
            Dists = cat(2,dist{:});
            ampRatio = max(AmpStarts(indx1),AmpEnds(indx2))./min(AmpStarts(indx1),AmpEnds(indx2));
            ampDiff = abs(AmpEnds(indx2)-AmpStarts(indx1));
            Skips = trackStartTime(indx2)-trackEndTime(indx1);
            
            %             costs = double(Dists.*(log2(ampRatio')+epsilon));
            %
            %condition on time
            IndsOfReasonableTimeJumps = logical((Skips>0).*(Skips<=maxTimeJump));
            indx1 = indx1(IndsOfReasonableTimeJumps);
            indx2 = indx2(IndsOfReasonableTimeJumps);
            Dists = Dists(IndsOfReasonableTimeJumps);
            ampRatio = ampRatio(IndsOfReasonableTimeJumps);
            ampDiff = ampDiff(IndsOfReasonableTimeJumps);
            Skips = Skips(IndsOfReasonableTimeJumps);
            clearvars IndsOfReasonableTimeJumps
            
            %condition on space
            IndsOfReasonableSpaceJumps = logical(Dists(:)<=sqrt(Skips(:))*maxStep);
            indx1 = indx1(IndsOfReasonableSpaceJumps);
            indx2 = indx2(IndsOfReasonableSpaceJumps);
            Dists = Dists(IndsOfReasonableSpaceJumps);
            ampRatio = ampRatio(IndsOfReasonableSpaceJumps);
            ampDiff = ampDiff(IndsOfReasonableSpaceJumps);
            Skips = Skips(IndsOfReasonableSpaceJumps);
            clearvars IndsOfReasonableSpaceJumps;
            
            
            %condition on intensity
            IndsOfReasonableIntJumps = logical(ampDiff(:)<=0.1);
            indx1 = indx1(IndsOfReasonableIntJumps);
            indx2 = indx2(IndsOfReasonableIntJumps);
            Dists = Dists(IndsOfReasonableIntJumps);
            ampRatio = ampRatio(IndsOfReasonableIntJumps);
            ampDiff = ampDiff(IndsOfReasonableIntJumps);
            Skips = Skips(IndsOfReasonableIntJumps);
            clearvars IndsOfReasonableIntJumps;
            
            
            %
            %             indx1=[];
            %             indx2=[];
            %             Dists=[];
            %             Skips = [];
            %
            %             f = waitbar(0,'Calculating ditances of possible links...');
            %             for ind1 = 1:numTracks
            %                 if ~rem(ind1,100)
            %                     waitbar(ind1/numTracks,f,'Calculating ditances of possible links...');
            %                 end
            %                 for ind2=1: numTracks
            %                     dT = trackStartTime(ind2)-trackEndTime(ind1);
            %                     if dT>0 && dT<=maxTimeJump;%condition on time
            %                         dR = norm(CentroidsStarts(ind1,:)-CentroidsEnds(ind2,:));
            %                         if dR<=sqrt(dT)*maxStep;%condition on space
            %                             indx1=[indx1 ind1];
            %                             indx2=[indx2 ind2];
            %                             Dists=[Dists dR];
            %                             Skips=[Skips dT];
            %                         end
            %                     end
            %                 end
            %             end
            %             close(f)
            %
            
            %% calculate cost matrix
            %f = waitbar(0,'Gap closing...');
            
            [trackStats,statsRelChange,errFlag] = getTrackStats(trackedFeatureInfo,0.70,maxTimeJump);
            dispSqTheta = trackStats.dispSqTheta;
            dispSqR = trackStats.dispSqR;
            addConst =   - dispSqR.*log(dispSqTheta) + log(gamma(dispSqR));%log(ampDiffStd)
            
            
            dispSqT = dispSqTheta(Skips);
            dispSqAR = dispSqR(Skips);
            addCons = addConst(Skips);
            costs = dispSqT(:).*Dists(:).^2-(dispSqAR(:)-1).*log(max(Dists(:).^2,realmin))+addCons(:);
            
            %             costs = [];
            %
            %             for i=1 : numel(indx1)
            %                 if ~rem(i,100)
            %                     waitbar(i/numel(indx1),f,'Calculating cost matrix for gap closing...');
            %                 end
            %                 dta = Skips(i);
            %                 distA = Dists(i);
            %                 dispSqT = dispSqTheta(dta);
            %                 dispSqAR = dispSqR(dta);
            %                 addCons = addConst(dta);
            %                 costs = [costs, dispSqT.*distA.^2-(dispSqAR-1).*log(max(distA.^2,realmin))+addCons];
            %             end
            %             close(f)
            
            
            
            
            costMat = sparse(indx1,indx2,costs,numTracks,numTracks);
            
            
            %% Close gaps with merging/splitting
            
            %if there are gaps to close (i.e. if there are tracks that start after the
            %first frame and tracks that end before the last frame) ...
            numTracksLink = size(trackedFeatureIndx,1);
            mergeSplit = 0;
            %if any(trackStartTime > 1) && any(trackEndTime < numFramesEff)
            
            
            
            %calculate the cost matrix, which already includes the
            %costs of birth and death
            % -- USER DEFINED FUNCTION -- %
            
            %link tracks based on this cost matrix, allowing for birth and death
            [link12,link21] = lap(costMat,[],[], 1, max(costs)+1);
            link12 = double(link12);
            link21 = double(link21);
            
            %put the indices of all tracks from linking in one vector
            tracks2Link = (1:numTracksLink)';
            tracksRemaining = tracks2Link;
            
            %reserve memory space for matrix showing track connectivity
            compoundTrack = zeros(numTracksLink,600);
            
            %initialize compTrackIndx
            compTrackIndx = 0;
            %
            while ~isempty(tracksRemaining)
                
                %update compound track index by 1
                compTrackIndx = compTrackIndx + 1
                
                %take first track as a seed to build a compound track with
                %closed gaps and merges/splits
                trackSeed = tracksRemaining(1);
                seedLength = 1;
                seedLengthOld = 0; %dummy just to get into the while loop
                
                %while current seed contains more tracks than previous seed, i.e.
                %whie new track segments are still being added to the compound
                %track
                while seedLength > seedLengthOld
                    
                    %store current seed for later comparison
                    seedLengthOld = seedLength;
                    
                    %find tracks connected to ends of seed tracks
                    tmpTracks = link12(trackSeed);
                    trackLink2End = tmpTracks(tmpTracks <= numTracksLink); %starts linked to ends
                    trackMerge = [];
                    if mergeSplit
                        trackMerge = indxMerge(tmpTracks(tmpTracks > numTracksLink & ...
                            tmpTracks <= numTracksLink+numMerge) - numTracksLink); %tracks that ends merge with
                    end
                    
                    %find tracks connected to starts of seed tracks
                    tmpTracks = link21(trackSeed);
                    trackLink2Start = tmpTracks(tmpTracks <= numTracksLink); %ends linked to starts
                    trackSplit = [];
                    if mergeSplit
                        trackSplit = indxSplit(tmpTracks(tmpTracks > numTracksLink & ...
                            tmpTracks <= numTracksLink+numSplit) - numTracksLink); %tracks that starts split from
                    end
                    
                    %put all tracks together as the new seed
                    trackSeed = [trackSeed; trackLink2End; trackLink2Start; ...
                        trackMerge; trackSplit];
                    
                    %remove repetitions and arrange tracks in ascending order
                    trackSeed = unique(trackSeed);
                    
                    %get number of tracks in new seed
                    seedLength = length(trackSeed);
                    
                    %expand new seed if merging/splitting are allowed
                    if mergeSplit
                        
                        %variables storing merge/split seed tracks
                        mergeSeed = [];
                        splitSeed = [];
                        
                        %go over all seed tracks
                        for iSeed = 1 : seedLength
                            
                            %get the location(s) of this track in indxMerge
                            mergeSeed = [mergeSeed; find(indxMerge == trackSeed(iSeed))];
                            
                            %get the location(s) of this track in indxSplit
                            splitSeed = [splitSeed; find(indxSplit == trackSeed(iSeed))];
                            
                        end
                        
                        %add numTracksLink to mergeSeed and splitSeed to determine
                        %their location in the cost matrix
                        mergeSeed = mergeSeed + numTracksLink;
                        splitSeed = splitSeed + numTracksLink;
                        
                        %find tracks merging with seed tracks
                        trackMerge = [];
                        for iSeed = 1 : length(mergeSeed)
                            trackMerge = [trackMerge; find(link12(1:numTracksLink)==mergeSeed(iSeed))];
                        end
                        
                        %find tracks splitting from seed tracks
                        trackSplit = [];
                        for iSeed = 1 : length(splitSeed)
                            trackSplit = [trackSplit; find(link21(1:numTracksLink)==splitSeed(iSeed))];
                        end
                        
                        %add these track to the seed
                        trackSeed = [trackSeed; trackMerge; trackSplit];
                        
                        %remove repetitions and arrange tracks in ascending order
                        trackSeed = unique(trackSeed);
                        
                        %get number of tracks in new seed
                        seedLength = length(trackSeed);
                        
                    end %(if mergeSplit)
                    
                end %(while length(trackSeed) > length(trackSeedOld))
                
                %expand trackSeed to reserve memory for connetivity information
                trackSeedConnect = [trackSeed zeros(seedLength,2)];
                
                %store the tracks that the ends of the seed tracks are linked to,
                %and indicate whether it's an end-to-start link (+ve) or a merge (-ve)
                tmpTracks = link12(trackSeed);
                if mergeSplit
                    tmpTracks(tmpTracks > numTracksLink & tmpTracks <= ...
                        numTracksLink+numMerge) = -indxMerge(tmpTracks(tmpTracks > ...
                        numTracksLink & tmpTracks <= numTracksLink+numMerge) - numTracksLink);
                end
                tmpTracks(tmpTracks > numTracksLink) = NaN;
                trackSeedConnect(:,2) = tmpTracks;
                
                %store the tracks that the starts of the seed tracks are linked to,
                %and indicate whether it's a start-to-end link (+ve) or a split (-ve)
                tmpTracks = link21(trackSeed);
                if mergeSplit
                    tmpTracks(tmpTracks > numTracksLink & tmpTracks <= ...
                        numTracksLink+numSplit) = -indxSplit(tmpTracks(tmpTracks > ...
                        numTracksLink & tmpTracks <= numTracksLink+numSplit) - numTracksLink);
                end
                tmpTracks(tmpTracks > numTracksLink) = NaN;
                trackSeedConnect(:,3) = tmpTracks;
                
                %store tracks making up this compound track and their connectivity
                compoundTrack(compTrackIndx,1:3*seedLength) = reshape(...
                    trackSeedConnect,3*seedLength,1)';
                
                %in the list of all tracks, indicate that these tracks have
                %been taken care of by placing NaN instead of their number
                tracks2Link(trackSeed) = NaN;
                
                %retain only tracks that have not been linked to anything yet
                tracksRemaining = tracks2Link(~isnan(tracks2Link));
                
            end %(while ~isempty(tracksRemaining))
            
            %remove empty rows
            maxValue = max(compoundTrack,[],2);
            compoundTrack = compoundTrack(maxValue > 0,:);
            
            %determine number of tracks after gap closing (including merge/split if
            %specified)
            numTracksCG = size(compoundTrack,1);
            
            %reserve memory for structure storing tracks after gap closing
            tracksFinal = repmat(struct('tracksFeatIndxCG',[],...
                'tracksCoordAmpCG',[],'seqOfEvents',[]),numTracksCG,1);
            
            f = waitbar(0,'Closing gaps...');
            
            %go over all compound tracks
            for iTrack = 1 : numTracksCG
                if ~rem(iTrack,100)
                    waitbar(iTrack/numTracksCG,f,'Closing gaps...');
                end
                %get indices of tracks from linking making up current compound track
                %determine their number and connectivity
                trackSeedConnect = compoundTrack(iTrack,:)';
                trackSeedConnect = trackSeedConnect(trackSeedConnect ~= 0);
                seedLength = length(trackSeedConnect)/3; %number of segments making current track
                trackSeedConnect = reshape(trackSeedConnect,seedLength,3);
                
                %get their start times
                segmentStartTime = trackStartTime(trackSeedConnect(:,1));
                
                %arrange segments in ascending order of their start times
                [segmentStartTime,indxOrder] = sort(segmentStartTime);
                trackSeedConnect = trackSeedConnect(indxOrder,:);
                
                %get the segments' end times
                segmentEndTime = trackEndTime(trackSeedConnect(:,1));
                
                %calculate the segments' positions in the matrix of coordinates and
                %amplitudes
                segmentStartTime8 = 8 * (segmentStartTime - 1) + 1;
                segmentEndTime8   = 8 * segmentEndTime;
                
                %instead of having the connectivity in terms of the original track
                %indices, have it in terms of the indices of this subset of tracks
                %(which are arranged in ascending order of their start times)
                for iSeed = 1 : seedLength
                    value = trackSeedConnect(iSeed,2);
                    if value > 0
                        trackSeedConnect(iSeed,2) = find(trackSeedConnect(:,1) == ...
                            value);
                    elseif value < 0
                        trackSeedConnect(iSeed,2) = -find(trackSeedConnect(:,1) == ...
                            -value);
                    end
                    value = trackSeedConnect(iSeed,3);
                    if value > 0
                        trackSeedConnect(iSeed,3) = find(trackSeedConnect(:,1) == ...
                            value);
                    elseif value < 0
                        trackSeedConnect(iSeed,3) = -find(trackSeedConnect(:,1) == ...
                            -value);
                    end
                end
                
                %get track information from the matrices storing linking information
                tracksFeatIndxCG = trackedFeatureIndx(trackSeedConnect(:,1),:);
                tracksCoordAmpCG = trackedFeatureInfo(trackSeedConnect(:,1),:);
                
                %convert zeros to NaNs where approriate for the case of sparse
                %matrices
                if issparse(tracksCoordAmpCG)
                    
                    %convert sparse to full
                    tracksCoordAmpCG = full(tracksCoordAmpCG);
                    
                    %go over all the rows in this compound track
                    for iRow = 1 : size(tracksCoordAmpCG,1)
                        
                        %find all the zero entries
                        colZero = find(tracksCoordAmpCG(iRow,:)==0);
                        colZero = colZero(:)';
                        
                        %find the columns of the x-coordinates corresponding to
                        %the zero columns
                        xCoordCol = colZero - mod(colZero-1,8*ones(size(colZero)));
                        
                        %keep only the columns whose x-coordinate is zero as
                        %well
                        colZero = colZero(tracksCoordAmpCG(iRow,xCoordCol)==0);
                        
                        %replace zero with NaN in the surviving columns
                        tracksCoordAmpCG(iRow,colZero) = NaN;
                        
                    end
                    
                end
                
                %perform all gap closing links and modify connectivity accordingly
                %go over all starts in reverse order
                for iSeed = seedLength : -1 : 2
                    
                    %find the track this track might be connected to
                    track2Append = trackSeedConnect(iSeed,3);
                    
                    %if there is a track (which is not a split)
                    if track2Append > 0
                        
                        %put track information in the relevant row
                        tracksFeatIndxCG(track2Append,segmentStartTime(iSeed):...
                            segmentEndTime(iSeed)) = tracksFeatIndxCG(iSeed,...
                            segmentStartTime(iSeed):segmentEndTime(iSeed));
                        tracksFeatIndxCG(iSeed,:) = 0;
                        tracksCoordAmpCG(track2Append,segmentStartTime8(iSeed):...
                            segmentEndTime8(iSeed)) = tracksCoordAmpCG(iSeed,...
                            segmentStartTime8(iSeed):segmentEndTime8(iSeed));
                        tracksCoordAmpCG(iSeed,:) = NaN;
                        
                        %update segment information
                        segmentEndTime(track2Append) = segmentEndTime(iSeed);
                        segmentEndTime8(track2Append) = segmentEndTime8(iSeed);
                        segmentEndTime(iSeed) = NaN;
                        segmentEndTime8(iSeed) = NaN;
                        segmentStartTime(iSeed) = NaN;
                        segmentStartTime8(iSeed) = NaN;
                        
                        %update connectivity
                        trackSeedConnect(track2Append,2) = trackSeedConnect(iSeed,2);
                        trackSeedConnect(trackSeedConnect(:,2) == iSeed,2) = track2Append;
                        trackSeedConnect(trackSeedConnect(:,3) == iSeed,3) = track2Append;
                        trackSeedConnect(trackSeedConnect(:,2) == -iSeed,2) = -track2Append;
                        trackSeedConnect(trackSeedConnect(:,3) == -iSeed,3) = -track2Append;
                        
                    end %(if track2Append > 0)
                    
                end %(for iSeed = seedLength : -1 : 2)
                
                %find rows that are not empty
                maxValue = max(tracksFeatIndxCG,[],2);
                rowsNotEmpty = find(maxValue > 0);
                
                %remove empty rows
                tracksFeatIndxCG = tracksFeatIndxCG(rowsNotEmpty,:);
                tracksCoordAmpCG = tracksCoordAmpCG(rowsNotEmpty,:);
                segmentEndTime   = segmentEndTime(rowsNotEmpty);
                segmentStartTime = segmentStartTime(rowsNotEmpty);
                trackSeedConnect = trackSeedConnect(rowsNotEmpty,:);
                
                %update connectivity accordingly
                %by now, only merges and splits are left - thus no need for minus
                %sign to distinguish them from closed gaps
                for iSeed = 1 : length(rowsNotEmpty)
                    trackSeedConnect(trackSeedConnect(:,2) == -rowsNotEmpty(...
                        iSeed),2) = iSeed;
                    trackSeedConnect(trackSeedConnect(:,3) == -rowsNotEmpty(...
                        iSeed),3) = iSeed;
                end
                
                %determine new "seedLength"
                seedLength = length(rowsNotEmpty);
                
                %store the sequence of events of this track
                seqOfEvents = [segmentStartTime ones(seedLength,1) ...
                    (1:seedLength)' trackSeedConnect(:,3); ...
                    segmentEndTime 2*ones(seedLength,1) ...
                    (1:seedLength)' trackSeedConnect(:,2)];
                
                %sort sequence of events in ascending order of time
                [tmp,indxOrder] = sort(seqOfEvents(:,1));
                seqOfEvents = seqOfEvents(indxOrder,:);
                
                %add 1 to the times of merges
                indx = find(~isnan(seqOfEvents(:,4)) & seqOfEvents(:,2) == 2);
                seqOfEvents(indx,1) = seqOfEvents(indx,1) + 1;
                
                %find the frame where the compound track starts and the frames
                %where it ends
                frameStart = seqOfEvents(1,1);
                frameEnd   = seqOfEvents(end,1);
                
                %store final tracks, removing frames before anything happens and
                %after everything happens
                tracksFinal(iTrack).tracksFeatIndxCG = tracksFeatIndxCG(:,...
                    frameStart:frameEnd);
                tracksFinal(iTrack).tracksCoordAmpCG = tracksCoordAmpCG(:,...
                    8*(frameStart-1)+1:8*frameEnd);
                tracksFinal(iTrack).seqOfEvents = seqOfEvents;
                
            end %(for iTrack = 1 : numTracksCG)
            close(f)
            %%
            
            a = arrayfun(@(x) size(x.tracksFeatIndxCG,2),tracksFinal);
            longtracksFinal = tracksFinal(a>15);% remove all tracks shorter than 15 frames
            
            R.Tracks = longtracksFinal;
        end
        
        
        
        function setEpitheliumFitParams(R,pos,varargin)
            CorneaCells = R.getCorneaCellsLbl(pos);
            
            for i=1:numel(CorneaCells)
                distScore = CorneaCells{i}.Centroids(:,3)-CorneaCells{i}.TopoZ;
                [h, Xbins] = histcounts(distScore,200,'Normalization', 'probability');
                Xbins = (Xbins(2:end)+Xbins(1:end-1))/2;
                if i==1 || nargin>2
                    plot(Xbins,h);
                    shg
                    title('select gaussian area for lower (epithelium) peak')
                    pause;
                    J =InAxes;
                else
                    J = logical((Xbins>(BETA(2)-2*BETA(3))).*(Xbins<(BETA(2)+BETA(3))));
                end
                hToFit = h(J);
                XtoFit = Xbins(J);
                
                %
                posit = median(XtoFit)
                stdev = std(XtoFit)
                amp = max(hToFit)*sqrt(2*pi)*stdev
                BETA0 = [amp posit stdev];
                [BETA,RESNORM,RESIDUAL,EXITFLAG] = lsqcurvefit(@GaussianFit, BETA0 ,XtoFit, hToFit,[0 -inf 0.8*stdev], [inf inf 1.2*stdev]);
                x = min(Xbins):0.1:max(Xbins);
                if EXITFLAG>=0;
                    plot(Xbins, h, '-.', x, GaussianFit(BETA, x));
                    figure(gcf)
                    if nargin>2
                        pause(0.2)
                    end
                end
                set(gca,'xlim',[-100,200],'ylim',[0,0.2]);
                drawnow;
                
                CorneaCells{i}.FitParam.func = @(x) GaussianFit(BETA,x);
                CorneaCells{i}.FitParam.cumFunc = @(x) (1+erf((x-BETA(2))/(sqrt(2)*BETA(3))))/2;
                CorneaCells{i}.FitParam.stdev = BETA(3);
                CorneaCells{i}.FitParam.posit = BETA(2);
                CorneaCells{i}.FitParam.amp = BETA(1);
                
                
                
                
                %Take upper and lower limits to keep as epithelium cells. +/-2\sigma
                lowLim = CorneaCells{i}.FitParam.posit-2*CorneaCells{i}.FitParam.stdev;
                highLim = CorneaCells{i}.FitParam.posit+2*CorneaCells{i}.FitParam.stdev;
                CorneaCells{i}.Jepi = ((distScore<highLim).*(distScore>lowLim));
                CorneaCells{i}.epiScore = 1-CorneaCells{i}.FitParam.cumFunc(distScore);
                CorneaCells{i}.epiScore(~CorneaCells{i}.Jepi)=NaN;
                
                CorneaCells{i}.Jepi = find(CorneaCells{i}.Jepi);
            end
        end
        
        function Data = loadRawStack(R,timepoint, tile, Channel, varargin)
            channelNames = R.partitionInfo.channelNames;
            indChNuc = strcmp(channelNames,Channel);
            partNameList = R.partitionInfo.name(find((R.partitionInfo.timepoint==timepoint).*(R.partitionInfo.tile==tile)));
            options = ParseInputs('options',  struct('level', 1,'waitbar',1), varargin);
            Data = loadBigDataViewerFormat([R.pth filesep partNameList{indChNuc}],options);
            Data = squeeze(single(Data)/2^12);
        end
        
        
        function R = merge(Rvec,varargin)
            % TODO: extend support for merging object with different
            % TimeVec structures, for now I'm assuming they are the same!
            arg.prefix='%g_';
            arg = parseVarargin(varargin,arg);
            arg.Results = MulticolorSPIMResults;
            R = merge@Results(Rvec,arg);
            
            R.SPIMTimepoints=SPIMTimepoints.empty(1,0);
            for i=1:numel(Rvec)
                R.PIVlbl = [R.SPIMTimepoints Rvec(i).SPIMTimepoints];
            end
            
            allfields = {};
            for i=1:numel(Rvec)
                allfields = union(allfields,fieldnames(Rvec(i).TimeVecs));
            end
            for i=1:numel(Rvec)
                T=Rvec(i).TimeVecs;
                missingflds = setdiff(allfields,fieldnames(T));
                for j=1:numel(missingflds)
                    T(1).(missingflds{j})=[];
                end
                T=orderfields(T,allfields);
                if i==1
                    TimeVec=T;
                else
                    TimeVec=[TimeVec T];  %#ok<AGROW>
                end
            end
            R.TimeVecs=TimeVec;
        end
        
        
        
        function h = plotCellVsTime(R,pos)
            CorneaCells = R.getCorneaCellsLbl(pos);
            plot((0:(numel(R.CorneaCellLbls{1})-1))/2,cellfun(@(x) x.num, CorneaCells));
            ylabel('#Cells')
            xlabel('Time(h)')
            shg
        end
        
        
        
        
        %plotting
        function HeatMapData(R, dataname,pos, frame, varargin)
            arg.clims = [-15, 15];
            arg = parseVarargin(varargin,arg);
            
            a = R.getData(dataname,pos);
            a = a(:,frame);
            if min(arg.clims)>=0
                colormap(magma())
            else
                colormap(makeColorMap([0.6 0 0.6],[0 0 0],[0.8 0.8 0]))
            end;
            imagesc(unique(R.PIVlbl{1}.X), unique(R.PIVlbl{1}.Y), reshape(a,numel(unique(R.PIVlbl{1}.X)),numel(unique(R.PIVlbl{1}.Y)))',arg.clims);
            set(gcf,'color','w');
            axis equal
            title(dataname)
            colorbar
            set(gca,'ydir','normal')
            shg;
        end
        
        
        function h = PlotDisp(R,pos, i, dt,varargin)
            
            %function to plot displacement bw frame i and i+dt
            CorneaCells = R.getCorneaCellsLbl(pos);
            CentroidsD = CorneaCells{i}.Centroids;
            CentroidsA = CorneaCells{i+dt}.Centroids;
            
            
            CM = double(CorneaCells{i}.Link12Mat);%continuous connectivity map
            for j=i+1:i+dt-1
                CM=CM*double(CorneaCells{j}.Link12Mat);
            end
            [indx1, indx2] = find(CM);
            
            [indx1, J] = sort(indx1);
            
            idxToPlot=J;
            
            
            indx2 = indx2(idxToPlot);
            CentroidsD = CentroidsD(indx1,:);
            CentroidsA = CentroidsA(indx2,:);
            
            range = 1:length(indx1);
            if nargin>4
                range = varargin{1};
                if range(end)>indx1(end)
                    range = range(1):indx1(end);
                    disp('range end out of bounds, drawing all points up to the total # of cells.')
                end
                if range(1)>indx1(end)
                    range = 1:length(indx1);
                    disp('range fully out of bounds, drawing all points.')
                end
            end;
            
            %tzeva = viridis(length(indx1));
            %scatter3(CentroidsD(range,1),CentroidsD(range,2),-CentroidsD(range,3),[],tzeva(range,:),'*');
            %hold on
            %scatter3(CentroidsA(range,1),CentroidsA(range,2),-CentroidsA(range,3),[],tzeva(range,:));
            dx=CentroidsA(range,1)-CentroidsD(range,1);
            dy=CentroidsA(range,2)-CentroidsD(range,2);
            dz=CentroidsA(range,3)-CentroidsD(range,3);
            
            h = quiver3(CentroidsD(range,1),CentroidsD(range,2),-CentroidsD(range,3),dx,dy,-dz, 0);
            
            
            set(gca,'xlim',[-200 3000],'ylim',[-200 3000],'CameraPositionMode','manual','CameraPosition',[-2.6364e+03 -1.4045e+04 2.1889e+03])
            %hold on
            shg
        end
        
        
        
        %% Some sanity checks
        function allTracks = allTrackMatrix(R)
            longtracksFinal = R.Tracks;
            frames = R.Frames;
            allTracks = zeros(size(longtracksFinal,1),numel(frames));
            for i=1:size(longtracksFinal,1);
                trackStart = longtracksFinal(i).seqOfEvents(1,1);
                trackEnd = longtracksFinal(i).seqOfEvents(end,1);
                allTracks(i,trackStart:trackEnd) = longtracksFinal(i).tracksFeatIndxCG;
            end
        end
        
        function allTracks = heatmapTracks(R,dataname)
            longtracksFinal = R.Tracks;
            frames = R.Frames;
            allTracks = zeros(size(longtracksFinal,1),numel(frames));
            for i=1:size(longtracksFinal,1);
                trackStart = longtracksFinal(i).seqOfEvents(1,1);
                trackEnd = longtracksFinal(i).seqOfEvents(end,1);
                allTracks(i,trackStart:trackEnd) = longtracksFinal(i).(dataname);
            end
        end
        
        function h = plotFracInTracks(R)
            h = plot(sum(allTrackMatrix(R)>0)./cellfun(@(x) x.num, R.SPIMTimepoints));
            xlabel('Frame')
            ylabel({'Fraction of detected cells' 'in an active track'})
            shg
        end
        
        function histTrackLength(R)
            histogram(sum(allTrackMatrix(R)>0,2));
            xlabel('Length (frames)')
            ylabel('#')
            shg
        end
        
        function idx = listTracksInRange(R, xRange,yRange,zRange)
            idx = arrayfun(@(x) R.findTrackInRange(x, xRange,yRange,zRange), R.Tracks);
        end
        
        function inOrNot = findTrackInRange(R,trk, xRange,yRange,zRange)
            inOrNot=1;
            inOrNot = inOrNot*(trk.tracksCoordAmpCG(1)>xRange(1));
            inOrNot = inOrNot*(trk.tracksCoordAmpCG(1)<xRange(2));
            inOrNot = inOrNot*(trk.tracksCoordAmpCG(2)>yRange(1));
            inOrNot = inOrNot*(trk.tracksCoordAmpCG(2)<yRange(2));
            inOrNot = inOrNot*(trk.tracksCoordAmpCG(3)>zRange(1));
            inOrNot = inOrNot*(trk.tracksCoordAmpCG(3)<zRange(2));
            inOrNot = logical(inOrNot);
        end
        
        
        function saveResults(R,pth)
            % save results to path pth
            if nargin==1
                pth = R.pth;
            end
            % make sure this is a pth based results (and not a merged
            % one...)
            assert(~isempty(pth),'can''t save Resutls with empty pth');
            
            R.SPIMTimepoints = {}; 
            
            % saving Conclusion as a seperate file to make it easier to read
            Conclusions = R.Conclusions;  %#ok<PROP,NASGU>
            

            save(fullfile(pth,'Results.mat'),'R','Conclusions', '-v7.3')
            fid=fopen(fullfile(pth,'Conclusions.txt'),'w');
            fprintf(fid,'%s',R.Conclusions);
            fclose(fid);
            save(fullfile(pth,'Conclusions.mat'),'Conclusions', '-v7.3')
            fid=fopen(fullfile(pth,'Classtype.txt'),'w');
            fprintf(fid,'%s',class(R));
            fclose(fid);
            if R.loadTPsflag
                SPIMTPs = cell(1,numel(R.Frames));
                for i=R.Frames'
                    disp(['loading timepoint ' num2str(i)])
                    SPIMTPs{i+1} = R.getSPIMTimepoints(i);
                end
                R.SPIMTimepoints = SPIMTPs;
            end
            
        end
        
    end
    
    
    
end