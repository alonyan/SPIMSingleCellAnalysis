function Blocklbl = SPIMBlockConstructor(fpath, partitionInfo,timepoint,tile)
    % set inits
    channelNames = partitionInfo.channelNames;
    NuclearChannel = 'DeepBlue';
    Blocklbl = SpimBlockLbl;
    Blocklbl.pth = fpath;
    Blocklbl.Frame = timepoint;
    Blocklbl.tile = tile;

    %find partitions that correspond to this position and timepoint
    partNameList = partitionInfo.name(find((partitionInfo.timepoint==timepoint).*(partitionInfo.tile==tile)));
    %make sure N files = N channels
    assert(numel(partNameList) == numel(partitionInfo.channelNames));
    Blocklbl.PartitionNameList = partNameList;
    
    
    %Load and process affine transforms
    tForms = partitionInfo.tForms(find((partitionInfo.timepoint==timepoint).*(partitionInfo.tile==tile).*(partitionInfo.channel==1)));   
    %Use associativity of affine transforms to create a single transform matrix; 
    T = eye(4);
    for i=1:numel(tForms{1})
        T = T*reshape([str2double(strsplit(tForms{1}{i}.affine.Text)), 0 0 0 1]',4,4); 
    end
    T = T';
    
    
    %Load channels
    options = struct('level', 1,'waitbar',0);
    %Channels = cell(numel(partNameList),1);
    Ints = cell(numel(partNameList),1);
    Ints90P = cell(numel(partNameList),1);
    
    %for i=1:numel(partNameList);
    %    Data = loadBigDataViewerFormat([fpath partNameList{i}],options);
    %    Data = squeeze(single(Data)/2^12);
    %    Channels{i} = Data;
    %end

    %Start calling cells, look at Nuc channel
    indChNuc = find(strcmp(channelNames,NuclearChannel));
    %Data = Channels{indChNuc};
    Data = loadBigDataViewerFormat([fpath partNameList{indChNuc}],options);
    Data = squeeze(single(Data)/2^12);
    Blocklbl.ImageDims = size(Data);

    % Gaussian smoothen and find regional max
    imgG = imgaussian3(Data,[5,1]);
    imgGRegMax = imregionalmax(imgG);
    % Expand regional max to a sphere of radius=4
    [xx,yy,zz] = ndgrid(-4:4);
    nhood = sqrt(xx.^2 + yy.^2+zz.^2) <= 4.0;
    ExpandedPeaks = imdilate(imgGRegMax,nhood);
    
    % Filter peaks with low signal using triangle thresholding
    somePix = Data(ExpandedPeaks);
    [h,x]=hist(log(somePix),100);
    ThreshVal = exp(triangleThresh(h,x));
    %ThreshVal = exp(meanThresh(log(somePix)));
    nuclearSeeds = ExpandedPeaks.*(imgG>ThreshVal);
    %nuclearSeedsSignal = Data.*nuclearSeeds;
    
    % Find connected components
    CC = bwconncomp(nuclearSeeds);
    
    
    % get region props in each channel
    S = regionprops(CC,Data,'MeanIntensity','Centroid','Area','PixelValues');
    % filter CC keep only cells with volume in some range
    Areas = cat(1, S.Area);
    Intensities = cat(1, S.MeanIntensity);
    Int90Prctile = arrayfun(@(x) prctile(x.PixelValues(x.PixelValues~=0),90),S);
    Centroids = cat(1,S.Centroid);
    %remove specks
    J = Areas>=120;
    Areas = Areas(J);
    Intensities = Intensities(J);
    Int90Prctile = Int90Prctile(J);
    Centroids = Centroids(J,:);
    
    %Add to Int cell array
    Ints{indChNuc}=Intensities;
    Ints90P{indChNuc}= Int90Prctile;

    %Now the other channels
    indOtherChannels = find(~strcmp(channelNames,NuclearChannel));
    for i=indOtherChannels
        %Data = Channels{i};
        Data = loadBigDataViewerFormat([fpath partNameList{indChNuc}],options);
        Data = squeeze(single(Data)/2^12);
        S = regionprops(CC,Data,'MeanIntensity','PixelValues');
        Intensities = cat(1, S.MeanIntensity);
        Int90Prctile = arrayfun(@(x) prctile(x.PixelValues(x.PixelValues~=0),90),S);
        %J is still the index of no specky cells
        Ints{i}=Intensities(J);
        Ints90P{i}= Int90Prctile(J);
    end
    
    %Apply transformations to centroids!
    tempCents = (T*[Centroids ones(1,size(Centroids,1))']')';
    Centroids = tempCents(:,1:3);
    
    Blocklbl.Centroids = Centroids;
    Blocklbl.Intensities = Ints;
    Blocklbl.Int90Prctile = Ints90P;

    Blocklbl.num = size(Centroids,1);
    %Blocklbl.DistFromManifold = DistFromManifold;
    %Blocklbl.TopoZ = TopoZ;
    %Blocklbl.GridZ = zz;
    Blocklbl.CC = CC;
    i
    
    
    
end
