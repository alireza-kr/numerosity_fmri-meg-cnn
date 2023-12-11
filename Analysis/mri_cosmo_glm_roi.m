%__________________________________________________________________________
function mri_cosmo_glm_roi(mypath,mri,mvpa,mode)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-ANALYSIS MODE - VOLUME
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(mode,'volume')
    
    warning('GLM is not implemented for Volume!');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-ANALYSIS MODE - SURFACE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(mode,'surface')

%-CHECK EXTERNALS
%==========================================================================
cosmo_check_external('surfing');
cosmo_check_external('afni');

%-IMPORT FILES AND INFORMATION
%==========================================================================
%-Target and Chunk
%targets = mvpa.cosmo.glm.targets;

%-Load predictor RDMs
if strcmp(mvpa.rdm.predictor.analysis,'Predictor')
    str.file = eval(mypath.file.data.rdm.predictor);
    str.folder = eval(mypath.folder.data.rdm);
    load(strcat(str.folder,'\',str.file));
    pred_cell = cell(2,length(mvpa.rdm.predictor.name));
    for i=1:length(mvpa.rdm.predictor.name)
       pred_cell{1,i} = mvpa.rdm.predictor.name{i};
       pred_cell{2,i} = eval(sprintf('RDM.%s.%s',...
           mvpa.rdm.predictor.type,...
           mvpa.rdm.predictor.name{i}));
    end

%-Load simulation RDMs
elseif strcmp(mvpa.rdm.predictor.analysis,'Simulation')
    str.file = eval(mypath.file.data.rdm.simulation);
    str.folder = eval(mypath.folder.data.rdm);
    load(strcat(str.folder,'\',str.file));
    pred_cell = cell(2,length(mvpa.rdm.predictor.name));
    for i=1:length(mvpa.rdm.predictor.name)
       pred_cell{1,i} = mvpa.rdm.predictor.name{i};
       pred_cell{2,i} = eval(sprintf('RDM.%s',...
           mvpa.rdm.predictor.name{i}));
    end
end

nsubjs = length(mri.subject.list);
nmasks = length(mvpa.cosmo.glm.mask.name);

%-Create output structure
table_size = [nsubjs nmasks];
table_type = repmat({'cell'},1,nmasks);
glm_table = table('Size',table_size,...
    'VariableTypes',table_type,...
    'VariableNames',mvpa.cosmo.glm.mask.name,...
    'RowNames',mri.subject.list);

%-LOAD DATA
%==========================================================================
data = load_data(mypath,mvpa,mri,'rdm');

%-ANALYSIS
%==========================================================================
for m = 1:nmasks
    for s = 1:nsubjs
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
        ds = data.(subj).ds.(mvpa.cosmo.glm.mask.name{m});
        %..................................................................
        
        %-ANALYSIS
        %------------------------------------------------------------------
        %-Regression
        measure = @MeasureRDMCorr;
        
        %-Measure arguments
        %..................................................................
        measure_args = struct();
        measure_args.mode = 'mri';
        measure_args.analysis = mvpa.cosmo.glm.analysis.mri;
        measure_args.input_dsm = true;
        measure_args.metric_dsm = mvpa.cosmo.glm.metric.mri;
        measure_args.type = mvpa.cosmo.glm.type.mri;
        measure_args.frrsa = mvpa.cosmo.glm.frrsa;
        measure_args.center_data = mvpa.cosmo.glm.voxscaling;
        measure_args.target_dsm = pred_cell(2,:);
        measure_args.labels = pred_cell(1,:)';
        %..................................................................
        
        %-Apply the measure
        ds_glm = measure(ds,measure_args);
        
        %-Beta values of each subject
        glm_table(s,m) = {num2str(ds_glm.samples)};
    end
end

%-PLOT THE RESULTS
%==========================================================================
for s = 1:nsubjs
    %-Create X and Y
    part = mvpa.cosmo.glm.mask.partition;
    x1 = categorical(mvpa.cosmo.glm.mask.name(1:part));
    x1 = reordercats(x1,mvpa.cosmo.glm.mask.name(1:part));
    x2 = categorical(mvpa.cosmo.glm.mask.name(part+1:end));
    x2 = reordercats(x2,mvpa.cosmo.glm.mask.name(part+1:end));
    for n = 1:length(ds_glm.sa.labels)
        for m = 1:length(mvpa.cosmo.glm.mask.name)
            sm_pred = str2num(glm_table{s,m}{:});
            y(n,m) = sm_pred(n);
        end
    end

    vis_args = struct();
    vis_args.part = part;
    vis_args.type = 'mri-subject-glm';
    vis_args.legend_name = pred_cell(1,:);
    vis_args.label_name.title = mvpa.vis.mri.glm.title;
    vis_args.label_name.xlabel = mvpa.vis.mri.glm.xlabel;
    vis_args.label_name.ylabel = mvpa.vis.mri.glm.ylabel;
    vis_args.color = cell2mat(values(mvpa.vis.color.glm,vis_args.legend_name)');
    DrawScatterSplit(x1,x2,y,vis_args);
    
    %-Save figure
    str.file = eval(mypath.file.figure.mri.glm.subj);
    str.folder = eval(mypath.folder.result.mri.figure);
    fn.figure.scatter = strcat(str.folder,'/',str.file);
    saveas(gcf,fn.figure.scatter);
    clf;
end

%-SAVE THE RESULTS
%==========================================================================
%-Save each subject result in a seperate mat file
for s = 1:nsubjs
    glm_table_subj = glm_table(s,:);
    str.file = eval(mypath.file.subj.mri.glm.roi);
    str.folder = eval(mypath.folder.result.mri.glm);
    fn.result.glm = strcat(str.folder,'/',str.file);
    save(fn.result.glm,'glm_table_subj','pred_cell');
end

close all;

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
function data = load_data(mypath,mvpa,mri,type)

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
