%__________________________________________________________________________
function noise_ceiling = MeasureNoiseCeiling(mypath,mri,meg,mvpa,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-This function uses RSA toolbox

%-Noise Ceiling Estimation
%https://www.johancarlin.com/understanding-noise-ceiling-metrics-rsa-compared-to-spearman-brown.html
%https://www.frontiersin.org/articles/10.3389/fncom.2015.00135/full
%https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003553

%-Studies used noise ceiling estimation for normalization
%https://direct.mit.edu/jocn/article/30/11/1559/28947/Tracking-the-Spatiotemporal-Neural-Dynamics-of
%https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008775
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

params=cosmo_structjoin('mode','mri',...
                        'sensor','',...
                        varargin);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-ANALYSIS MODE - MRI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(params.mode,'mri')

%-CHECK EXTERNALS
%==========================================================================
cosmo_check_external('surfing');
cosmo_check_external('afni');

%-IMPORT FILES AND INFORMATION
%==========================================================================
%-Target and Chunk
targets = mvpa.cosmo.glm.targets;

nsubjs = length(mri.subject.list);
nmasks = length(mvpa.cosmo.glm.mask.name);

table_size = [2 nmasks];
table_type = repmat({'double'},1,nmasks);
noise_ceiling = table('Size',table_size,...
    'VariableTypes',table_type,...
    'VariableNames',mvpa.cosmo.glm.mask.name);

%-LOAD DATA
%==========================================================================
data = load_data_mri(mypath,mvpa,mri,'rdm');

%-ANALYSIS
%==========================================================================
for m=1:nmasks
    refRDMestimates = {};
    
    for s=1:nsubjs
        subj = mri.subject.list{s};
        
        %-LOAD SUBJECT DATA
        %------------------------------------------------------------------
        %-Raw
        %..................................................................
        %ds_full = read_surf(mypath,mri,mvpa,subj,m,data);
        %ds_full.sa.targets = targets;
        %-Remove constant features
        %ds = cosmo_remove_useless_data(ds_full);
        %..................................................................

        %-RDM
        %..................................................................
        ds_dsm = data.(subj).ds.(mvpa.cosmo.glm.mask.name{m});
        %..................................................................
        
        %-COMPUTE RDM
        %--------------------------------------------------------------------
        %..................................................................
        %params.metric = mvpa.cosmo.glm.metric.mri;
        %params.center_data = mvpa.cosmo.glm.voxscaling;
        %params.pseudo = 0;
        %params.perm = 0;
        %params.cv = 'none';
        %params.mnn = 'none';
        %..................................................................
        %mode = strcat(params.metric,'-',params.mnn,'-',params.cv,'-','mri');
        %ds_dsm = MeasureNeuralRDM(ds,mode,params);
        
        refRDMestimates(s).name = subj;
        refRDMestimates(s).RDM = ds_dsm;
        
    end
    
    [ceiling_upperBound,ceiling_lowerBound,~] = ceilingAvgRDMcorr(refRDMestimates,mvpa.cosmo.glm.type.mri);
    noise_ceiling(1,m) = {ceiling_upperBound};
    noise_ceiling(2,m) = {ceiling_lowerBound};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-ANALYSIS MODE - MEG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(params.mode,'meg')

%-LOAD DATA
%==========================================================================
data = load_data_meg(mypath,mvpa,meg,params.sensor,'rdm');

%-ANALYSIS
%==========================================================================
nsubjs = length(meg.subject.list);
for s = 1:nsubjs
    subj = meg.subject.list{s};

    %-LOAD SUBJECT DATA
    %----------------------------------------------------------------------
    %-Raw
    %......................................................................
    %ds = data.(subj).ds;
    %......................................................................

    %-RDM
    %......................................................................
    ds_dsm = data.(subj).ds;
    %......................................................................

    %-COMPUTE RDM
    %--------------------------------------------------------------------
    %....................................................................
    %params.metric = mvpa.cosmo.glm.metric.meg;
    %params.center_data = mvpa.cosmo.glm.chanscaling;
    %params.pseudo = mvpa.cosmo.glm.npseudo;
    %params.perm = mvpa.cosmo.glm.nperm;
    %params.cv = mvpa.cosmo.glm.cv;
    %params.mnn = mvpa.cosmo.glm.mnn;
    %....................................................................
    %mode = strcat(params.metric,'-',params.mnn,'-',params.cv,'-','meg');
    %ds_dsm = MeasureNeuralRDM(ds,mode,params);

    for t=1:meg.epoch.ntimes
        refRDMestimates(s).name = subj;
        refRDMestimates(s).RDM = ds_dsm(:,:,t);
        refRDM{t} = refRDMestimates;
    end
end

for t=1:meg.epoch.ntimes
    [ceiling_upperBound,ceiling_lowerBound,~] = ceilingAvgRDMcorr(refRDM{t},mvpa.cosmo.glm.type.meg);
    noise_ceiling(1,t) = {ceiling_upperBound};
    noise_ceiling(2,t) = {ceiling_lowerBound};
end

end

end

%__________________________________________________________________________
function ds_full = read_surf(mypath,mri,mvpa,subj,m,data)

    ds_full = struct();
    
    activeVoxels = read_voxel(mypath,mri,mvpa,subj,m,data);
    
    for s = mri.analysis.split.num
        for d=1:length(mvpa.cosmo.glm.file.(mri.task.name))
            %-Select active voxels
            ds = cosmo_slice(data.(subj){s,d},logical(activeVoxels),2);

            if(isequal(ds_full,struct()))
                ds_full=ds;
            else
                ds_full = cosmo_stack({ds_full,ds});
            end
        end
    end
    
end

%__________________________________________________________________________
function activeVoxels = read_voxel(mypath,mri,mvpa,subj,m,data)

    %-If the desired number of voxels within each ROIs is set to be 0,
    %then select all possible voxels within the ROI.
    if mri.mask.nvox==0
        %-Read the atlas file
        activeVoxels = data.mask{m};
        activeVoxels = activeVoxels(:,2)';
    else
        %-Read the atlas file
        mask = data.mask{m};
        
        %-Read the source file
        source = data.source{find(strcmp(mri.subject.list,subj))};
        
        %-Looking for voxels that are in the mask and have a nonzero t-value
        inmask=find(source.samples~=0 & mask(:,2)'>0);
        mask_vol_nonzeroVox=zeros(1,length(source.samples));
        mask_vol_nonzeroVox(inmask)=source.samples(inmask);
        [sorted,inds]=sort(mask_vol_nonzeroVox(inmask));
        if length(sorted) < mri.mask.nvox
            disp('Insufficient number of voxels in mask!');
        end
        mask_vol_activeVox=zeros(1,length(source.samples));
        mask_vol_activeVox(inmask(inds(length(sorted)-(mri.mask.nvox-1):length(sorted))))=1;
        source.samples = mask_vol_activeVox;
            
        activeVoxels = source.samples;
        
    end
end

%__________________________________________________________________________
function data = load_data_mri(mypath,mvpa,mri,type)

nsubjs = length(mri.subject.list);
nmasks = length(mvpa.cosmo.glm.mask.name);

if strcmp(type,'raw')
    %-Read Atlas
    %----------------------------------------------------------------------
    str.file = eval(mypath.file.data.atlas);
    str.folder = eval(mypath.folder.data.atlas);
    fn.atlas = strcat(str.folder,'/',str.file);
    [err,msk,Info,Com] = Read_1D(fn.atlas);

    for m = 1:length(mvpa.cosmo.glm.mask.name) 
        data.mask{m}(:,1) = msk(:,1);
        %-Create the mask by putting to one the index of the ROIs of interest
        %and putting to zero the values of other regions.
        data.mask{m}(:,2) = ismember(msk(:,2),mvpa.cosmo.glm.mask.map(mvpa.cosmo.glm.mask.name{m}));
    end

    %-Read Sources
    %----------------------------------------------------------------------
    for s = 1:nsubjs
        subj = mri.subject.list{s};

        str.file = eval(mypath.file.subj.mri.source.(mri.mask.voxel.type).dset);
        str.folder = eval(mypath.folder.result.mri.source.(mri.mask.voxel.type));
        fn.source = strcat(str.folder,'/',str.file);
        if strcmp(mri.mask.voxel.type,'firstLevel')
            data.source{s} = cosmo_surface_dataset(fn.source);
        elseif strcmp(mri.mask.voxel.type,'secondLevel')
            data.source{s} = cosmo_surface_dataset(fn.source);
            data.source{s}.samples = data.source{s}.samples(1,:);
        end
    end

    %-Read Subjects' Data (1st-Level T-Map)
    %----------------------------------------------------------------------
    for sbj = 1:nsubjs
        subj = mri.subject.list{sbj};
        for s = mri.analysis.split.num
            for d=1:length(mvpa.cosmo.glm.file.(mri.task.name))
                c = mvpa.cosmo.glm.file.(mri.task.name){d};
                str.file = eval(mypath.file.subj.mri.firstLevel.dset);
                str.folder = eval(mypath.folder.subj.mri.firstLevel);
                fn.data = strcat(str.folder,'/',str.file);
                data.(subj){s,d} = cosmo_surface_dataset(fn.data);
            end
        end
    end
    
elseif strcmp(type,'rdm')
    for s = 1:nsubjs
        subj = mri.subject.list{s};
        
        str.file = eval(mypath.file.subj.mri.rdm);
        str.folder = eval(mypath.folder.result.mri.rdm);
        fn.in = strcat(str.folder,'/',str.file);
        load(fn.in);
        data.(subj).ds = RDM;
    end
end

end

%__________________________________________________________________________
function data = load_data_meg(mypath,mvpa,meg,sensor,type)

nsubjs = length(meg.subject.list);

if strcmp(type,'raw')
    for s = 1:nsubjs
        subj = meg.subject.list{s};
        
        str.file = eval(mypath.file.subj.meg.preprocessed);
        str.folder = eval(mypath.folder.subj.meg);
        fn.in = strcat(str.folder,'/',str.file);
        data.(subj).ds = load(fn.in);

        %-Load data of Gradiometer or Magnetometer based on sensor parameter.
        if strcmp(sensor,'mag')
            data.(subj).ds.trial = data.(subj).ds.trial(:,data.(subj).ds.channel_type(:,1)=='m',:);
        elseif strcmp(sensor,'grad')
            data.(subj).ds.trial = data.(subj).ds.trial(:,data.(subj).ds.channel_type(:,1)=='g',:);
        else
            warning("Please choose sensor type (mag or grad)!");
        end
    end
elseif strcmp(type,'rdm')
    for s = 1:nsubjs
        subj = meg.subject.list{s};
        
        str.file = eval(mypath.file.subj.meg.rdm);
        str.folder = eval(mypath.folder.result.meg.rdm);
        fn.in = strcat(str.folder,'/',str.file);
        load(fn.in);
        data.(subj).ds = RDM.(sensor);
    end
end

end
