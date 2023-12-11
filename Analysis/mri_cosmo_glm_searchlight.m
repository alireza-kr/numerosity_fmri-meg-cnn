%http://cosmomvpa.org/ex_surface_searchlight.html
%http://www.cosmomvpa.org/_static/publish/run_rsm_measure_searchlight.html
%https://www.cosmomvpa.org/_static/publish/run_surface_searchlight.html

%__________________________________________________________________________
function mri_cosmo_glm_searchlight(mypath,mri,mvpa,mode)

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
targets = mvpa.cosmo.glm.targets;

%-Load predictor RDMs
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

nsubjs = length(mri.subject.list);

%-ANALYSIS
%==========================================================================
for s = 1:nsubjs
    subj = mri.subject.list{s};

    %-IMPORT SURFACES
    %======================================================================
    str.file = eval(mypath.file.surface.pial);
    str.folder = eval(mypath.folder.subj.surface.surfing);
    fn.pial = strcat(str.folder,'/',str.file);
    str.file = eval(mypath.file.surface.wm);
    fn.white = strcat(str.folder,'/',str.file);

    [pial_v,pial_f] = surfing_read(fn.pial);
    fprintf('The pial surface has %d vertices, %d faces\n',size(pial_v,1),size(pial_f,1))

    [white_v,white_f] = surfing_read(fn.white);
    fprintf('The white surface has %d vertices, %d faces\n',size(pial_v,1),size(pial_f,1))

    %-verify that the face information in pial_f and white_f are the same
    assert(isequal(white_f,pial_f));

    %-show the content of the surfaces
    fprintf('pial_v\n');
    cosmo_disp(pial_v)
    fprintf('pial_f\n');
    cosmo_disp(pial_f)
    fprintf('white_v\n');
    cosmo_disp(white_v)
    fprintf('white_f\n');
    cosmo_disp(white_f)

    %-COMPUTE THICKNESS OF THE CORTEX
    %======================================================================
    delta = pial_v - white_v;
    delta_squared = delta .^ 2;
    thickness_squared = sum(delta_squared,2);
    thickness = sqrt(thickness_squared);

    %-READ INFALTED SURFACE
    %======================================================================
    %-For visualization purposes, read inflated surface
    str.file = eval(mypath.file.surface.inflated);
    str.folder = eval(mypath.folder.subj.surface.surfing);
    fn.inflated = strcat(str.folder,'/',str.file);
    [infl_v,infl_f]=surfing_read(fn.inflated);
    fprintf('The inflated surface has %d vertices, %d faces\n',size(infl_v,1),size(infl_f,1))
    
    nvertices=size(infl_v,1);

    %-SAVE THICKNESS SURFACE
    %======================================================================
    ds_thickness=struct();
    ds_thickness.fa.node_indices=1:nvertices;
    ds_thickness.samples=thickness(:)';
    ds_thickness.a.fdim.labels={'node_indices'};
    ds_thickness.a.fdim.values={(1:nvertices)'};

    str.file = eval(mypath.file.surface.thickness);
    str.folder = eval(mypath.folder.result.mri.searchlight);
    fn.result.thickness=strcat(str.folder,'/',str.file);
    cosmo_map2surface(ds_thickness,fn.result.thickness);
    
    %-LOAD DATA SET
    %======================================================================
    ds_full = read_surf(mypath,mri,mvpa,subj);
    ds_full.sa.targets = targets;
    
    %-ANALYSIS
    %======================================================================
    %-Do not remove constant features as you get the following error:
    %https://groups.google.com/g/cosmomvpa/c/FiPjOsA3nmI/m/b3YxOe5nBQAJ?pli=1
    ds = cosmo_remove_useless_data(ds_full);
    
    %-Regression
    measure = @MeasureRDMCorr;

    %-Measure arguments
    %......................................................................
    measure_args = struct();
    measure_args.mode = 'mri';
    measure_args.analysis = mvpa.cosmo.glm.analysis.mri;
    measure_args.metric_dsm = mvpa.cosmo.glm.metric.mri;
    measure_args.type = mvpa.cosmo.glm.type.mri;
    measure_args.frrsa = mvpa.cosmo.glm.frrsa;
    measure_args.center_data = mvpa.cosmo.glm.voxscaling;
    measure_args.target_dsm = pred_cell(2,:);
    measure_args.labels = pred_cell(1,:)';
    %......................................................................
    
    surface_def={pial_v,white_v,pial_f};

    %nbrhood=cosmo_surficial_neighborhood(ds,surface_def,'count',mvpa.cosmo.searchlight.nvox);
    nbrhood=cosmo_surficial_neighborhood(ds,surface_def,'radius',mvpa.cosmo.searchlight.r);
    
    %-Run searchlight
    ds_sl=cosmo_searchlight(ds,nbrhood,measure,measure_args);
    
    %-SAVE FILE
    %======================================================================
    for p=1:length(pred_cell(1,:))
        str.file = eval(mypath.file.subj.mri.glm.searchlight);
        str.folder = eval(mypath.folder.result.mri.searchlight);
        fn.result.glm = strcat(str.folder,'/',str.file);
        cosmo_map2surface(cosmo_slice(ds_sl,p),fn.result.glm);
    end
    
end

end

end

%__________________________________________________________________________
function ds_full = read_surf(mypath,mri,mvpa,subj)

    ds_full = struct();

    for s = mri.analysis.split.num
        for d=1:length(mvpa.cosmo.glm.file.(mri.task.name))
            c = mvpa.cosmo.glm.file.(mri.task.name){d};
            str.file = eval(mypath.file.subj.mri.firstLevel.dset);
            str.folder = eval(mypath.folder.subj.mri.firstLevel);
            fn.data = strcat(str.folder,'/',str.file);
            data = cosmo_surface_dataset(fn.data);

            if(isequal(ds_full,struct()))
                ds_full = data;
            else
                ds_full = cosmo_stack({ds_full,data});
            end
        end
    end
    
end
