%% Load MD
BaseStr = '/RazorScopeData/RazorScopeSets/';
Usr = 'Zach';
Project = 'CorneaHSV';
Dataset = 'Infection72H_2019May21';
acquisition = 1;
%% Get MD of raw data
acqname = ['acq_' num2str(acquisition)];
fpath = [BaseStr Usr filesep Project filesep Dataset filesep acqname ];
%% Get metadata

MD=Metadata(fpath);

channelNames = arrayfun(@(x) MD.getSpecificMetadataByIndex('Channel', x), 1:numel(MD.unique('Channel')));
info = xml2struct([fpath filesep 'hdf5_dataset.xml']);
partitionInfo = struct;
partitionInfo.('id') = cellfun(@(x) str2double(x.setups.Text),info.SpimData.SequenceDescription.ImageLoader.partition);
partitionInfo.('name') = cellfun(@(x) x.path.Text,info.SpimData.SequenceDescription.ImageLoader.partition, 'UniformOutput',false);
partitionInfo.('timepoint') = cellfun(@(x) str2double(x.timepoints.Text),info.SpimData.SequenceDescription.ImageLoader.partition);
partitionInfo.('tForms') = cellfun(@(x) x.ViewTransform, info.SpimData.ViewRegistrations.ViewRegistration, 'UniformOutput',false);


beadInfo = struct;
beadInfo.('filename') = cellfun(@(x) x.Text, info.SpimData.ViewInterestPoints.ViewInterestPointsFile, 'UniformOutput',false);
beadInfo.('id') = cellfun(@(x) x.Attributes.setup, info.SpimData.ViewInterestPoints.ViewInterestPointsFile, 'UniformOutput',false);
beadInfo.('timepoint') = cellfun(@(x) x.Attributes.timepoint, info.SpimData.ViewInterestPoints.ViewInterestPointsFile, 'UniformOutput',false);
beadInfo.('label') = cellfun(@(x) x.Attributes.label, info.SpimData.ViewInterestPoints.ViewInterestPointsFile, 'UniformOutput',false);
%make sure to only take the relevant interest points
J = find(cellfun(@(x) strcmp(x,'beads'), beadInfo.label));
beadInfo.filename = {beadInfo.filename{J}};
beadInfo.id = {beadInfo.id{J}};
beadInfo.timepoint = {beadInfo.timepoint{J}};
beadInfo.label = {beadInfo.label{J}};

partitionInfo.beadsFile = cell(size(partitionInfo.id));

ids = str2double(beadInfo.id);
tps = str2double(beadInfo.timepoint);

for i=unique(tps)
    for j=unique(ids)
        partitionInfo.beadsFile{logical((partitionInfo.id==j).*(partitionInfo.timepoint==i))}= [beadInfo.filename{logical((ids==j).*(tps==i))} '.ip.txt'];
    end
end




setupInfo.id = cellfun(@(x) str2double(x.id.Text), info.SpimData.SequenceDescription.ViewSetups.ViewSetup);
setupInfo.tile = cellfun(@(x) str2double(x.attributes.tile.Text), info.SpimData.SequenceDescription.ViewSetups.ViewSetup);
setupInfo.channel = cellfun(@(x) str2double(x.attributes.channel.Text), info.SpimData.SequenceDescription.ViewSetups.ViewSetup);
setupInfo.angle = cellfun(@(x) str2double(x.attributes.angle.Text), info.SpimData.SequenceDescription.ViewSetups.ViewSetup);

partitionInfo.tile = zeros(size(partitionInfo.id));
partitionInfo.channel = zeros(size(partitionInfo.id));
partitionInfo.angle = zeros(size(partitionInfo.id));
for i=setupInfo.id
    partitionInfo.tile(partitionInfo.id==i) = setupInfo.tile(setupInfo.id==i);
    partitionInfo.channel(partitionInfo.id==i) = setupInfo.channel(setupInfo.id==i);
    partitionInfo.angle(partitionInfo.id==i) = setupInfo.angle(setupInfo.id==i);
end
clear setupInfo
partitionInfo.channelNames = channelNames;



%% Make results object

R = MulticolorSPIMResults(fpath)
R.partitionInfo = partitionInfo;
%frames = unique(cell2mat(MD.getSpecificMetadata('frame')));
R.Frames = 0:72;
R.analysisScript=fullfile([fpath filesep 'AnalysisScriptTemplate.m']);%Change this to the right file
R.reportPth = [BaseStr 'Reports' filesep 'Alon' filesep Project filesep Dataset];
R.NuclearChannel = 'Red';

% define nuclear channel and channel to use for tracking
NucChannel = 'Red'; %Channel to segment on
TrackChannel = 'Red'; %"smoothness" channel, usually virus
%R.saveResults

%% Big calculation. Find single cells, etc.
for i=[0:72]
    TP = R.getSPIMTimepoints(i);
end


%SPIMTimepoints = cell(numel(R.Frames),1);
%%

%% Load MD
BaseStr = '/RazorScopeData/RazorScopeSets/';
Usr = 'Zach';
Project = 'CorneaHSV';
Dataset = 'Infection72H_2019May21';
acquisition = 1;
%% Get MD of raw data
acqname = ['acq_' num2str(acquisition)];
fpath = [BaseStr Usr filesep Project filesep Dataset filesep acqname ];
%% Get metadata
MD=Metadata(fpath);
R = MulticolorSPIMResults(fpath)


%% Link adjecent frames
TrackChannel = 'Red'
R.Link(TrackChannel)


%% Close gaps
R.closeGaps(TrackChannel)
R.saveResults








%% Calculate more features about single cell tracks


for WellNum=1:numel(Wells)
    Tracks = R.getTracks(R.PosNames{WellNum});
    
    %Time
    T = arrayfun(@(x) (x.seqOfEvents(1,1):x.seqOfEvents(2,1)), Tracks,'UniformOutput',false);
    [Tracks.('T')] = T{:};
    
    
    WellCells = R.getWellsLbl(R.PosNames{WellNum});
    
    indtrckChnl = find(strcmp(TrackChannel,WellCells{1}.channels));
    indNucChnl = find(strcmp(NucChannel,WellCells{1}.channels));
    
    
    VirusIntensities = cellfun(@(x) x.Int90Prctile{indtrckChnl}, WellCells,'UniformOutput', false);
    NuclearIntensities = cellfun(@(x) x.Int90Prctile{indNucChnl}, WellCells,'UniformOutput', false);
    
    for i=1:numel(Tracks)
        i
        VirusTrack{i} = zeros(1,numel(Tracks(i).T));
        NuclearTrack{i} = zeros(1,numel(Tracks(i).T));
        
        for j=1:numel(Tracks(i).T)
            % j
            if Tracks(i).tracksFeatIndxCG(j)
                VirusTrack{i}(j) = VirusIntensities{Tracks(i).T(j)}(Tracks(i).tracksFeatIndxCG(j))';
                NuclearTrack{i}(j) = NuclearIntensities{Tracks(i).T(j)}(Tracks(i).tracksFeatIndxCG(j))';
            else
                VirusTrack{i}(j) = NaN;
                NuclearTrack{i}(j) = NaN;
            end
        end
    end
    [Tracks.('VirusTrack')] = VirusTrack{:};
    [Tracks.('NuclearTrack')] = NuclearTrack{:};
    
    
    WellCells = R.getWellsLbl(R.PosNames{WellNum});
    ThreshVirus = mean(WellCells{1}.Int90Prctile{indtrckChnl})+2*std(WellCells{1}.Int90Prctile{indtrckChnl});
    ThreshNuc = mean(WellCells{1}.Int90Prctile{indNucChnl})+2*std(WellCells{1}.Int90Prctile{indNucChnl});
    
    for i=1:numel(Tracks)
        Tracks(i).Infected = Smoothing(Tracks(i).VirusTrack)>ThreshVirus;
        Tracks(i).CellsGetInfected = sum(Tracks(i).Infected)>4;
        Tracks(i).Dead = Smoothing(Tracks(i).NuclearTrack,'neigh', 5)>ThreshNuc;
        Tracks(i).CellDies = sum(Tracks(i).Dead)>8  && any(diff(Tracks(i).Dead)==1);% && sum(diff(Tracks(i).Dead))==1;
    end
    
    R.setTracks(Tracks,R.PosNames{WellNum})
end
clearvars NuclearTrack  VirusTrack Tracks;
R.saveResults

