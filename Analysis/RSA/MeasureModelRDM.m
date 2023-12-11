%__________________________________________________________________________
function MeasureModelRDM(mypath,mri,meg,mvpa,type)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%http://www.cosmomvpa.org/matlab/cosmo_dissimilarity_matrix_measure.html
%http://www.cosmomvpa.org/_static/publish/run_rsa_visualize.html

%-DISTATIS
%https://www.hindawi.com/journals/cmmm/2013/796183/
%https://www.cosmomvpa.org/matlab/demo_fmri_distatis.html

%-Create RDM with SVM
%https://groups.google.com/g/cosmomvpa/c/QfWg-VcJdXg
%https://osf.io/u67z3/
%https://osf.io/2637j/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

param = rdm(type);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-Predictor RDM (Castaldi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(param.rdm.type,'Castaldi')

%..........................................................................
diff_str_1 = 'RDM.Diff.%s(i,j) = abs(log10(param.rdm.stim.%s(i)) - log10(param.rdm.stim.%s(j)));';
diff_str_2 = 'RDM.Diff.%s(j,i) = RDM.Diff.%s(i,j);';

%-Create Matrixes
for dim=1:length(param.rdm.dimension.list)
    for i=1:param.rdm.size
        for j=i:param.rdm.size
            %-Distance: Difference
            eval(sprintf(diff_str_1,...
                param.rdm.dimension.list{dim},...
                param.rdm.dimension.list{dim},...
                param.rdm.dimension.list{dim}));
            eval(sprintf(diff_str_2,...
                param.rdm.dimension.list{dim},...
                param.rdm.dimension.list{dim}));
        end
    end
end

%..........................................................................
norm_diff_str_1 = 'RDM.NormDiff.%s = RDM.Diff.%s - min(RDM.Diff.%s(:));';
norm_diff_str_2 = 'RDM.NormDiff.%s = RDM.NormDiff.%s ./ max(RDM.NormDiff.%s(:));';

%-Normalize Matrixes between 0 to 1
for dim=1:length(param.rdm.dimension.list)
    eval(sprintf(norm_diff_str_1,...
        param.rdm.dimension.list{dim},...
        param.rdm.dimension.list{dim},...
        param.rdm.dimension.list{dim}));
    eval(sprintf(norm_diff_str_2,...
        param.rdm.dimension.list{dim},...
        param.rdm.dimension.list{dim},...
        param.rdm.dimension.list{dim}));
end

%..........................................................................
for i=1:param.rdm.size
    for j=1:length(param.rdm.dimension.list)
        %-i: Stimulus (32), j: Features (N,S,TFA,TSA,Dens)
        RDM.features(i,j) = param.rdm.stim.(param.rdm.dimension.list{j})(i);
    end
end

%-ZCA-Whitening the RDMs
%--------------------------------------------------------------------------
fields = fieldnames(RDM.Diff);
for f=1:numel(fields)
    vec_rdm(:,f) = cosmo_squareform(RDM.Diff.(fields{f}),'tovector')';
end

[vec_rdm, ~, ~, ~] = whiten(vec_rdm,0.0001);

for f=1:numel(fields)
    RDM.WDiff.(fields{f}) = cosmo_squareform(vec_rdm(:,f));
end
%..........................................................................
fields = fieldnames(RDM.NormDiff);
for f=1:numel(fields)
    vec_rdm(:,f) = cosmo_squareform(RDM.NormDiff.(fields{f}),'tovector')';
end

[vec_rdm, ~, ~, ~] = whiten(vec_rdm,0.0001);

for f=1:numel(fields)
    RDM.WNormDiff.(fields{f}) = cosmo_squareform(vec_rdm(:,f));
end

str.file = eval(mypath.file.data.rdm.predictor);
str.folder = eval(mypath.folder.data.rdm);
fn.result.rdm = strcat(str.folder,'\',str.file);
save(fn.result.rdm,'RDM');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-Predictor RDM (DeWind)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(param.rdm.type,'DeWind')

%..........................................................................
diff_str_1 = 'RDM.Diff.%s(i,j) = abs(log10(param.rdm.stim.%s(i)) - log10(param.rdm.stim.%s(j)));';
diff_str_2 = 'RDM.Diff.%s(j,i) = RDM.Diff.%s(i,j);';

%-Create Matrixes
for dim=1:length(param.rdm.dimension.list)
    for i=1:param.rdm.size
        for j=i:param.rdm.size
            %-Distance: Difference
            eval(sprintf(diff_str_1,...
                param.rdm.dimension.list{dim},...
                param.rdm.dimension.list{dim},...
                param.rdm.dimension.list{dim}));
            eval(sprintf(diff_str_2,...
                param.rdm.dimension.list{dim},...
                param.rdm.dimension.list{dim}));
        end
    end
end

%..........................................................................
norm_diff_str_1 = 'RDM.NormDiff.%s = RDM.Diff.%s - min(RDM.Diff.%s(:));';
norm_diff_str_2 = 'RDM.NormDiff.%s = RDM.NormDiff.%s ./ max(RDM.NormDiff.%s(:));';

%-Normalize Matrixes between 0 to 1
for dim=1:length(param.rdm.dimension.list)
    eval(sprintf(norm_diff_str_1,...
        param.rdm.dimension.list{dim},...
        param.rdm.dimension.list{dim},...
        param.rdm.dimension.list{dim}));
    eval(sprintf(norm_diff_str_2,...
        param.rdm.dimension.list{dim},...
        param.rdm.dimension.list{dim},...
        param.rdm.dimension.list{dim}));
end

str.file = eval(mypath.file.data.rdm.predictor);
str.folder = eval(mypath.folder.data.rdm);
fn.result.rdm = strcat(str.folder,'\',str.file);
save(fn.result.rdm,'RDM');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-Simulation RDM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(param.rdm.type,'Simulation')
    
ntargets = length(mvpa.sim.num.condition); nchunks = 1;
ds_size = mvpa.sim.ds.size;

%-Monotonic (Square Root)
%......................................................................
%-Dataset
for n=1:length(mvpa.sim.num.unique)
    ds_num{n} = cosmo_synthetic_dataset('size',ds_size,'ntargets',ntargets,'nchunks',nchunks);
    for t=1:ntargets
        ds_num{n}.samples(t,:) = sqrt(mvpa.sim.num.condition(t));
    end
end
ds{1} = cosmo_stack(ds_num,2);

%-RDM
params.metric = mvpa.cosmo.glm.metric.mri;
params.center_data = mvpa.cosmo.glm.voxscaling;
params.pseudo = 0;
params.perm = 0;
params.cv = 'none';
params.mnn = 'none';
mode = strcat(params.metric,'-',params.mnn,'-',params.cv,'-','mri');
RDM.(mvpa.sim.func.name{1}) = MeasureNeuralRDM(ds{1},mode,params);
%......................................................................

%-Monotonic (Logarithm)
%......................................................................
%-Dataset
for n=1:length(mvpa.sim.num.unique)
    ds_num{n} = cosmo_synthetic_dataset('size',ds_size,'ntargets',ntargets,'nchunks',nchunks);
    for t=1:ntargets
        ds_num{n}.samples(t,:) = log(mvpa.sim.num.condition(t));
    end
end
ds{2} = cosmo_stack(ds_num,2);

%-RDM
params.metric = mvpa.cosmo.glm.metric.mri;
params.center_data = mvpa.cosmo.glm.voxscaling;
params.pseudo = 0;
params.perm = 0;
params.cv = 'none';
params.mnn = 'none';
mode = strcat(params.metric,'-',params.mnn,'-',params.cv,'-','mri');
RDM.(mvpa.sim.func.name{2}) = MeasureNeuralRDM(ds{2},mode,params);
%......................................................................

%-Tunning 1 (SD: 2)
%......................................................................
%-Dataset
for n=1:length(mvpa.sim.num.unique)
    ds_num{n} = cosmo_synthetic_dataset('size',ds_size,'ntargets',ntargets,'nchunks',nchunks);
    sd = mvpa.sim.func.sd(3);
    mean_val = mvpa.sim.num.unique(n);
    for t=1:ntargets
        x = mvpa.sim.num.condition(t);
        ds_num{n}.samples(t,:) = 1/(2*pi*sd)*exp(-(x-mean_val).^2/(2*sd^2));
    end
end
ds{3} = cosmo_stack(ds_num,2);

%-RDM
params.metric = mvpa.cosmo.glm.metric.mri;
params.center_data = mvpa.cosmo.glm.voxscaling;
params.pseudo = 0;
params.perm = 0;
params.cv = 'none';
params.mnn = 'none';
mode = strcat(params.metric,'-',params.mnn,'-',params.cv,'-','mri');
RDM.(mvpa.sim.func.name{3}) = MeasureNeuralRDM(ds{3},mode,params);
%......................................................................

%-Tunning 2 (SD: 4)
%......................................................................
%-Dataset
for n=1:length(mvpa.sim.num.unique)
    ds_num{n} = cosmo_synthetic_dataset('size',ds_size,'ntargets',ntargets,'nchunks',nchunks);
    sd = mvpa.sim.func.sd(4);
    mean_val = mvpa.sim.num.unique(n);
    for t=1:ntargets
        x = mvpa.sim.num.condition(t);
        ds_num{n}.samples(t,:) = 1/(2*pi*sd)*exp(-(x-mean_val).^2/(2*sd^2));
    end
end
ds{4} = cosmo_stack(ds_num,2);

%-RDM
params.metric = mvpa.cosmo.glm.metric.mri;
params.center_data = mvpa.cosmo.glm.voxscaling;
params.pseudo = 0;
params.perm = 0;
params.cv = 'none';
params.mnn = 'none';
mode = strcat(params.metric,'-',params.mnn,'-',params.cv,'-','mri');
RDM.(mvpa.sim.func.name{4}) = MeasureNeuralRDM(ds{4},mode,params);
%......................................................................

%-Tunning 3 (SD: 6)
%......................................................................
%-Dataset
for n=1:length(mvpa.sim.num.unique)
    ds_num{n} = cosmo_synthetic_dataset('size',ds_size,'ntargets',ntargets,'nchunks',nchunks);
    sd = mvpa.sim.func.sd(5);
    mean_val = mvpa.sim.num.unique(n);
    for t=1:ntargets
        x = mvpa.sim.num.condition(t);
        ds_num{n}.samples(t,:) = 1/(2*pi*sd)*exp(-(x-mean_val).^2/(2*sd^2));
    end
end
ds{5} = cosmo_stack(ds_num,2);

%-RDM
params.metric = mvpa.cosmo.glm.metric.mri;
params.center_data = mvpa.cosmo.glm.voxscaling;
params.pseudo = 0;
params.perm = 0;
params.cv = 'none';
params.mnn = 'none';
mode = strcat(params.metric,'-',params.mnn,'-',params.cv,'-','mri');
RDM.(mvpa.sim.func.name{5}) = MeasureNeuralRDM(ds{5},mode,params);
%......................................................................

str.file = eval(mypath.file.data.rdm.simulation);
str.folder = eval(mypath.folder.data.rdm);
fn.result.rdm = strcat(str.folder,'\',str.file);
save(fn.result.rdm,'RDM');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-GIST RDM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%https://people.csail.mit.edu/torralba/code/spatialenvelope/
elseif strcmp(param.rdm.type,'GIST')
    
    warning('GIST RDM is not implemented yet!');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-DNN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%https://www.frontiersin.org/articles/10.3389/fninf.2021.679838/full
%https://github.com/ViCCo-Group/THINGSvision
elseif strcmp(param.rdm.type,strcat(mvpa.dnn.architecture,'-',mvpa.dnn.database))

for l=1:length(mvpa.dnn.layer.name)
%-Load data set
%--------------------------------------------------------------------------
str.file = eval(mypath.file.data.dnn);
str.folder = eval(mypath.folder.data.dnn.feature);
fn.in = strcat(str.folder,'/',str.file);
layer_python = single(readNPY(fn.in));

%-Average instances
%--------------------------------------------------------------------------
for t=1:length(mvpa.dnn.input.instant.type)
    layer(t,:) = mean(layer_python(mvpa.dnn.input.instant.num*(t-1)+1:mvpa.dnn.input.instant.num*t,:),1);
end
clear layer_python;

%-Analysis
%--------------------------------------------------------------------------

%......................................................................
params.metric = mvpa.cosmo.glm.metric.dnn;
params.center_data = false;
params.cv = 'none';
params.mnn = 'none';
%......................................................................
mode = strcat(params.metric,'-',params.mnn,'-',params.cv,'-','dnn');
ds_dsm = MeasureNeuralRDM(layer,mode,params);

RDM.(matlab.lang.makeValidName(mvpa.dnn.layer.name{l})) = ds_dsm;

%-Visualize RDM
%......................................................................
vis_args = struct();
vis_args.type =mvpa.dnn.layer.name{l};
vis_args.label = mvpa.dnn.input.instant.type;
DrawRDM(ds_dsm,vis_args);

str.file = eval(mypath.file.figure.dnn.rdm);
str.folder = eval(mypath.folder.result.dnn.figure);
fn.figure.rdm = strcat(str.folder,'/',str.file);
saveas(gcf,fn.figure.rdm);
clf;

%-Visualize MDS
%......................................................................
vis_args = struct();
vis_args.cof_num = mvpa.cosmo.decoding.cof.num;
vis_args.cof_inst = mvpa.cosmo.decoding.cof.inst;
vis_args.label = mvpa.cosmo.decoding.label;
vis_args.color = distinguishable_colors(mvpa.cosmo.decoding.cof.num);
DrawMDS(ds_dsm,vis_args);

str.file = eval(mypath.file.figure.dnn.mds);
str.folder = eval(mypath.folder.result.dnn.figure);
fn.figure.mds = strcat(str.folder,'/',str.file);
saveas(gcf,fn.figure.mds);
clf;

clear layer;

end

%-Save RDM
%--------------------------------------------------------------------------
str.file = eval(mypath.file.data.rdm.dnn);
str.folder = eval(mypath.folder.data.rdm);
fn.result.rdm = strcat(str.folder,'/',str.file);
save(fn.result.rdm,'RDM');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-Null DNN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(param.rdm.type,strcat('null-',mvpa.dnn.architecture,'-',mvpa.dnn.database))

for l=1:length(mvpa.dnn.layer.name)
%-Load data set
%--------------------------------------------------------------------------
str.file = eval(mypath.file.data.dnn);
str.folder = eval(mypath.folder.data.dnn.feature);
fn.in = strcat(str.folder,'/',str.file);
layer_python = single(readNPY(fn.in));

%-Average instances
%--------------------------------------------------------------------------
for t=1:length(mvpa.dnn.input.instant.type)
    layer(t,:) = mean(layer_python(mvpa.dnn.input.instant.num*(t-1)+1:mvpa.dnn.input.instant.num*t,:),1);
end
clear layer_python;

%-Create Null Distribution
%--------------------------------------------------------------------------
permuted_order_mat = [];
for sz=1:mvpa.dnn.null.size

    %-Shuffle Rows (Labels of data)
    %----------------------------------------------------------------------
    permutation_not_accepted = true;
    while permutation_not_accepted
        permuted_order = randperm(size(layer, 1));
        isequal_n = isequal(mvpa.rdm.type.N,mvpa.rdm.type.N(permuted_order));
        isequal_s = isequal(mvpa.rdm.type.S,mvpa.rdm.type.S(permuted_order));
        isequal_tfa = isequal(mvpa.rdm.type.TFA,mvpa.rdm.type.TFA(permuted_order));
        isequal_tsa = isequal(mvpa.rdm.type.TSA,mvpa.rdm.type.TSA(permuted_order));
        isequal_dens = isequal(mvpa.rdm.type.Dens,mvpa.rdm.type.Dens(permuted_order));
        if (~isequal_n) && (~isequal_s) && (~isequal_tfa) && (~isequal_tsa) && (~isequal_dens)
            if isequal(permuted_order_mat,[])
                permuted_order_mat = [permuted_order_mat;permuted_order];
                permutation_not_accepted = false;
            elseif ~ismember(permuted_order,permuted_order_mat,'rows')
                permuted_order_mat = [permuted_order_mat;permuted_order];
                permutation_not_accepted = false;
            end
        end
    end
    
end

for sz=1:mvpa.dnn.null.size
    
    permuted_layer = layer(permuted_order_mat(sz,:),:);

    %-Analysis
    %----------------------------------------------------------------------

    %......................................................................
    params.metric = mvpa.cosmo.glm.metric.dnn;
    params.center_data = false;
    params.cv = 'none';
    params.mnn = 'none';
    %......................................................................
    mode = strcat(params.metric,'-',params.mnn,'-',params.cv,'-','dnn');
    ds_dsm = MeasureNeuralRDM(permuted_layer,mode,params);

    NULL{sz}.RDM.(matlab.lang.makeValidName(mvpa.dnn.layer.name{l})) = ds_dsm;

end

clear layer;

end

%-Save RDM
%--------------------------------------------------------------------------
for sz=1:mvpa.dnn.null.size
    str.file = eval(mypath.file.data.null.dnn);
    str.folder = eval(mypath.folder.result.dnn.null);
    fn.result.rdm = strcat(str.folder,'/',str.file);
    RDM = NULL{sz}.RDM;
    save(fn.result.rdm,'RDM');
    clear RDM
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-Feature Reweighted DNN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%https://www.frontiersin.org/articles/10.3389/fninf.2021.679838/full
%https://github.com/ViCCo-Group/THINGSvision
elseif strcmp(param.rdm.type,strcat('fr-',mvpa.dnn.architecture,'-',mvpa.dnn.database))

%-Load Group-Average MRI RDM
%--------------------------------------------------------------------------
str.file.mri = eval(mypath.file.data.rdm.avg.mri);
str.folder = eval(mypath.folder.data.rdm);
load(strcat(str.folder,'\',str.file.mri),'RDM');
RDM_MRI = RDM;
clear RDM;

for l=1:length(mvpa.dnn.layer.name)
%-Load data set
%--------------------------------------------------------------------------
str.file = eval(mypath.file.data.dnn);
str.folder = eval(mypath.folder.data.dnn.feature);
fn.in = strcat(str.folder,'/',str.file);
layer_python = single(readNPY(fn.in));

%-Average instances
%--------------------------------------------------------------------------
for t=1:length(mvpa.dnn.input.instant.type)
    layer(t,:) = mean(layer_python(mvpa.dnn.input.instant.num*(t-1)+1:mvpa.dnn.input.instant.num*t,:),1);
end
clear layer_python;

%-Analysis
%--------------------------------------------------------------------------
lyr = matlab.lang.makeValidName(mvpa.dnn.layer.name{l});

for m = 1:length(mvpa.cosmo.glm.mask.name)
    msk = mvpa.cosmo.glm.mask.name{m};

    %-COMPUTE RDM
    %------------------------------------------------------------------
    %-Define targer
    target = RDM_MRI.(msk);
    %-Define predictor
    predictor = layer;

    RDM.(msk).(lyr) = MeasureFRRSA(mypath,target,predictor,'rdm');

    %-Visualize RDM
    %......................................................................
    vis_args = struct();
    vis_args.type = mvpa.dnn.layer.name{l};
    vis_args.label = mvpa.dnn.input.instant.type;
    DrawRDM(RDM.(msk).(lyr),vis_args);

    %str.file = eval();
    str.folder = eval(mypath.folder.result.dnn.figure);
    fn.figure.rdm = strcat(str.folder,'/',str.file);
    saveas(gcf,fn.figure.rdm);
    clf;

    %-Visualize MDS
    %......................................................................
    vis_args = struct();
    vis_args.cof_num = mvpa.cosmo.decoding.cof.num;
    vis_args.cof_inst = mvpa.cosmo.decoding.cof.inst;
    vis_args.label = mvpa.cosmo.decoding.label;
    vis_args.color = distinguishable_colors(mvpa.cosmo.decoding.cof.num);
    DrawMDS(RDM.(msk).(lyr),vis_args);

    %str.file = eval();
    str.folder = eval(mypath.folder.result.dnn.figure);
    fn.figure.mds = strcat(str.folder,'/',str.file);
    saveas(gcf,fn.figure.mds);
    clf;
end

clear layer;

end

%-Save RDM
%--------------------------------------------------------------------------
%str.file = eval();
str.folder = eval(mypath.folder.data.rdm);
fn.result.rdm = strcat(str.folder,'/',str.file);
save(fn.result.rdm,'RDM');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-Group-Average MRI RDM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(param.rdm.type,'MRI')

%-CHECK EXTERNALS
%--------------------------------------------------------------------------
cosmo_check_external('surfing');
cosmo_check_external('afni');
    
%-IMPORT FILES AND INFORMATION
%--------------------------------------------------------------------------
%-Target and Chunk
targets = mvpa.cosmo.glm.targets;

nsubjs = length(mri.subject.list);
nmasks = length(mvpa.cosmo.glm.mask.name);

%-LOAD DATA
%==========================================================================
data = load_data_mri(mypath,mvpa,mri);

%-Load predictor RDMs to be regressed out
reg_out_rdm_name = {'TSA','Dens'};
str.file = eval(mypath.file.data.rdm.predictor);
str.folder = eval(mypath.folder.data.rdm);
load(strcat(str.folder,'\',str.file));
for i=1:length(reg_out_rdm_name)
   reg_out_rdm = eval(sprintf('RDM.%s.%s',...
       mvpa.rdm.predictor.type,...
       reg_out_rdm_name{i}));
   reg_out_vec(:,i) = cosmo_squareform(reg_out_rdm,'tovector')';
end

%-ANALYSIS
%==========================================================================
for m = 1:nmasks
    
    %-Sum RDM
    rdm_sum = zeros();
    
    for s = 1:nsubjs
        subj = mri.subject.list{s};
        
        %-Load data set
        ds_full = read_surf(mypath,mri,mvpa,subj,m,data);
        ds_full.sa.targets = targets;
        
        %-COMPUTE RDM
        %------------------------------------------------------------------
        %-Compute average for each unique target, so that the datasets
        %have one sample for each target
        ds = cosmo_fx(ds_full,@(x)mean(x,1),'targets',1);
        
        %-Remove constant features
        ds = cosmo_remove_useless_data(ds);
        
        %..................................................................
        params.metric = mvpa.cosmo.glm.metric.mri;
        params.center_data = mvpa.cosmo.glm.voxscaling;
        params.pseudo = 0;
        params.perm = 0;
        params.cv = 'none';
        params.mnn = 'none';
        %..................................................................
        mode = strcat(params.metric,'-',params.mnn,'-',params.cv,'-','mri');
        ds_dsm = MeasureNeuralRDM(ds,mode,params);
        
        Subject_RDM.(subj).(mvpa.cosmo.glm.mask.name{m}) = ds_dsm;
        rdm_sum = rdm_sum + ds_dsm;
    end
    
    %-Average RDM
    Average_RDM.(mvpa.cosmo.glm.mask.name{m}) = rdm_sum/nsubjs;
    
    %-Visualize RDM
    vis_args = struct();
    vis_args.type = mvpa.cosmo.glm.mask.name{m};
    vis_args.label = mvpa.rdm.data.name;
    DrawRDM(Average_RDM.(mvpa.cosmo.glm.mask.name{m}),vis_args);
    
    str.file = eval(mypath.file.figure.mri.rdm.avg);
    str.folder = eval(mypath.folder.result.mri.figure);
    fn.figure.rdm = strcat(str.folder,'/',str.file);
    saveas(gcf,fn.figure.rdm);
    clf;
    
    %-Visualize Dendrogram
    
    %-Visualize MDS
    vis_args = struct();
    vis_args.cof_num = mvpa.cosmo.decoding.cof.num;
    vis_args.cof_inst = mvpa.cosmo.decoding.cof.inst;
    vis_args.label = mvpa.cosmo.decoding.label;
    vis_args.color = distinguishable_colors(mvpa.cosmo.decoding.cof.num);
    DrawMDS(Average_RDM.(mvpa.cosmo.glm.mask.name{m}),vis_args);
    
    str.file = eval(mypath.file.figure.mri.mds.avg);
    str.folder = eval(mypath.folder.result.mri.figure);
    fn.figure.mds = strcat(str.folder,'/',str.file);
    saveas(gcf,fn.figure.mds);
    clf;
    
    %-Average Number RDM
    Average_Num_RDM.(mvpa.cosmo.glm.mask.name{m}) = ...
        makeNumRDM(mvpa,Average_RDM.(mvpa.cosmo.glm.mask.name{m}),reg_out_vec,false);
    
    %-Visualize Number RDM
    vis_args = struct();
    vis_args.type = mvpa.cosmo.glm.mask.name{m};
    vis_args.label = mvpa.cosmo.tf.num.label;
    DrawRDM(Average_Num_RDM.(mvpa.cosmo.glm.mask.name{m}),vis_args);
    
    str.file = eval(mypath.file.figure.mri.rdm_num.avg);
    str.folder = eval(mypath.folder.result.mri.figure);
    fn.figure.rdm_num = strcat(str.folder,'/',str.file);
    saveas(gcf,fn.figure.rdm_num);
    clf;
    
    %-Tuning Function
    %[indiv,merged] = makeTF(mvpa,Average_Num_RDM.(mvpa.cosmo.glm.mask.name{m}));
    
    %-Visualize Gaussian Function of Numbers (Individual)
    %vis_args = struct();
    %vis_args.func = 'individual';
    %vis_args.type = mvpa.cosmo.glm.mask.name{m};
    %vis_args.label = mvpa.cosmo.tf.num.label;
    %vis_args.color = distinguishable_colors(mvpa.cosmo.tf.cof.num);
    %DrawFunction(indiv.f,indiv.x,indiv.y,indiv.fwhm,vis_args);
    
    %str.file = eval(mypath.file.figure.mri.tf.avg.indiv);
    %str.folder = eval(mypath.folder.result.mri.figure);
    %fn.figure.tf = strcat(str.folder,'/',str.file);
    %saveas(gcf,fn.figure.tf);
    %clf;
    
    %-Visualize Gaussian Function of Numbers (Merged)
    %vis_args = struct();
    %vis_args.func = 'merged';
    %vis_args.type = mvpa.cosmo.glm.mask.name{m};
    %DrawFunction(merged.f,merged.x,merged.y,merged.fwhm,vis_args);
    
    %str.file = eval(mypath.file.figure.mri.tf.avg.merged);
    %str.folder = eval(mypath.folder.result.mri.figure);
    %fn.figure.tf = strcat(str.folder,'/',str.file);
    %saveas(gcf,fn.figure.tf);
    %clf;
    
end

%-Save Subject RDM
%..........................................................................
for s = 1:nsubjs
    subj = mri.subject.list{s};
    
    str.file = eval(mypath.file.subj.mri.rdm);
    str.folder = eval(mypath.folder.result.mri.rdm);
    fn.result.rdm = strcat(str.folder,'/',str.file);
    RDM = Subject_RDM.(subj);
    save(fn.result.rdm,'RDM');
end
%..........................................................................

%-Save Average RDM
%..........................................................................
str.file = eval(mypath.file.data.rdm.avg.mri);
str.folder = eval(mypath.folder.data.rdm);
fn.result.rdm = strcat(str.folder,'/',str.file);
RDM = Average_RDM;
RDM.info.subjects = mri.subject.list;
save(fn.result.rdm,'RDM');
%..........................................................................

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-Group-Average MRI RDM (ROI)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(param.rdm.type,'MRI-ROI')

%-CHECK EXTERNALS
%--------------------------------------------------------------------------
cosmo_check_external('surfing');
cosmo_check_external('afni');

%-LOAD DATA
%==========================================================================
mri.mask.atlas = 'wang-ventral';
str.file = eval(mypath.file.data.rdm.avg.mri);
str.folder = eval(mypath.folder.data.rdm);
fn.result.rdm = strcat(str.folder,'/',str.file);
ventral = load(fn.result.rdm,'RDM');

mri.mask.atlas = 'wang-dorsal';
str.file = eval(mypath.file.data.rdm.avg.mri);
str.folder = eval(mypath.folder.data.rdm);
fn.result.rdm = strcat(str.folder,'/',str.file);
dorsal = load(fn.result.rdm,'RDM');

rois.evc = {'V1','V2','V3'};
rois.dorsal = {'V3AB','V7','IPS12','IPS345'};
rois.ventral = {'hV4','VO1','VO2','PHC1','PHC2'};
rois.all = {'V1','V2','V3','V3AB','V7','IPS12','IPS345','hV4','VO1','VO2','PHC1','PHC2'};

for r=1:length(rois.evc)
    RDM.(rois.evc{r}) = dorsal.RDM.(rois.evc{r});
end

for r=1:length(rois.dorsal)
    RDM.(rois.dorsal{r}) = dorsal.RDM.(rois.dorsal{r});
end

for r=1:length(rois.ventral)
    RDM.(rois.ventral{r}) = ventral.RDM.(rois.ventral{r});
end

%-CREATE DATASET
%==========================================================================
for r=1:length(rois.all)
    ds.samples(r,:) = cosmo_squareform(RDM.(rois.all{r}),'tovector');
end

ds.sa.targets = 1:length(rois.all);

%-ANALYSIS
%==========================================================================
%..................................................................
params.metric = mvpa.cosmo.glm.metric.mri;
params.center_data = mvpa.cosmo.glm.voxscaling;
params.pseudo = 0;
params.perm = 0;
params.cv = 'none';
params.mnn = 'none';
%..................................................................
mode = strcat(params.metric,'-',params.mnn,'-',params.cv,'-','mri');
ds_dsm = MeasureNeuralRDM(ds,mode,params);

%-VISUALIZE & SAVE RDM
%==========================================================================
%-Visualize RDM
%..........................................................................
vis_args = struct();
vis_args.type = 'ROIs';
vis_args.label = rois.all;
DrawRDM(ds_dsm,vis_args);

str.file = eval(mypath.file.figure.mri.rdm_roi.avg);
str.folder = eval(mypath.folder.result.mri.figure);
fn.figure.rdm = strcat(str.folder,'/',str.file);
saveas(gcf,fn.figure.rdm);

%-Save RDM
%..........................................................................
str.file = eval(mypath.file.data.rdm.avg.roi);
str.folder = eval(mypath.folder.data.rdm);
fn.result.rdm = strcat(str.folder,'/',str.file);
RDM = ds_dsm;
save(fn.result.rdm,'RDM');
%..........................................................................

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-Group-Average MEG RDM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(param.rdm.type,'MEG')

%-IMPORT FILES AND INFORMATION
%--------------------------------------------------------------------------
nsubjs = length(meg.subject.list);
sensors.type = {'grad','mag'};
sensors.name = {'Gradiometer','Magnetometer'};

%-ANALYSIS
%--------------------------------------------------------------------------
%-Iterate through sensors
for sen=1:length(sensors.type)
    sensor = sensors.type{sen};
    
    %-Sum RDM
    rdm_sum = zeros();
    
    for s = 1:nsubjs
        subj = meg.subject.list{s};

        %-Load data set
        str.file = eval(mypath.file.subj.meg.preprocessed);
        str.folder = eval(mypath.folder.subj.meg);
        fn.in = strcat(str.folder,'/',str.file);
        ds = load(fn.in);

        %-Load data of Gradiometer or Magnetometer based on
        %sensor parameter.
        if strcmp(sensor,'mag')
            ds.trial = ds.trial(:,ds.channel_type(:,1)=='m',:);
        elseif strcmp(sensor,'grad')
            ds.trial = ds.trial(:,ds.channel_type(:,1)=='g',:);
        else
            warning("Please choose sensor type (mag or grad)!");
        end

        %-COMPUTE RDM
        %------------------------------------------------------------------
        %..................................................................
        params.metric = mvpa.cosmo.glm.metric.meg;
        params.center_data = mvpa.cosmo.glm.chanscaling;
        params.pseudo = mvpa.cosmo.glm.npseudo;
        params.perm = mvpa.cosmo.glm.nperm;
        params.cv = mvpa.cosmo.glm.cv;
        params.mnn = mvpa.cosmo.glm.mnn;
        %..................................................................
        mode = strcat(params.metric,'-',params.mnn,'-',params.cv,'-','meg');
        ds_dsm = MeasureNeuralRDM(ds,mode,params);
        
        Subject_RDM.(subj).(sensor) = ds_dsm;
        rdm_sum = rdm_sum + ds_dsm;
    end
    
    %-Average RDM
    Average_RDM.(sensor) = rdm_sum/nsubjs;
end

%-Save Subject RDM
%..........................................................................
for s = 1:nsubjs
    subj = meg.subject.list{s};
    
    str.file = eval(mypath.file.subj.meg.rdm);
    str.folder = eval(mypath.folder.result.meg.rdm);
    fn.result.rdm = strcat(str.folder,'/',str.file);
    RDM = Subject_RDM.(subj);
    save(fn.result.rdm,'RDM');
end
%..........................................................................

%-Save Average RDM
%..........................................................................
str.file = eval(mypath.file.data.rdm.avg.meg);
str.folder = eval(mypath.folder.data.rdm);
fn.result.rdm = strcat(str.folder,'/',str.file);
RDM = Average_RDM;
RDM.info.subjects = meg.subject.list;
save(fn.result.rdm,'RDM');
%..........................................................................

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-Group-Average MEG-Frequency RDM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(param.rdm.type,'MEG-Frequency')

%-IMPORT FILES AND INFORMATION
%--------------------------------------------------------------------------
nsubjs = length(meg.subject.list);
sensors.type = {'grad','mag'};
sensors.name = {'Gradiometer','Magnetometer'};

%-ANALYSIS
%--------------------------------------------------------------------------
%-Iterate through sensors
for sen=1:length(sensors.type)
    sensor = sensors.type{sen};
    
    %-Sum RDM
    rdm_sum = cell(1,meg.convo.frex_num);
    rdm_sum(1,:) = {zeros()};
    
    for s = 1:nsubjs
        subj = meg.subject.list{s};
        
        %-Load data set
        str.file = eval(mypath.file.subj.meg.preprocessed);
        str.folder = eval(mypath.folder.subj.meg);
        fn.in = strcat(str.folder,'/',str.file);
        ds=load(fn.in);
        
        %-Load data of Gradiometer or Magnetometer based on
        %sensor parameter.
        if strcmp(sensor,'mag')
            ds.trial = ds.trial(:,ds.channel_type(:,1)=='m',:);
        elseif strcmp(sensor,'grad')
            ds.trial = ds.trial(:,ds.channel_type(:,1)=='g',:);
        else
            warning("Please choose sensor type (mag or grad)!");
        end
        
        %-DECOMPOSE TO FREQUENCY DOMAIN
        %------------------------------------------------------------------
        [ds.freq, ds.convo] = timefrexdecomp(meg,ds);
        
        for f=1:meg.convo.frex_num
            
            dsFreq.trial = cat(1,ds.freq{:,f});
            dsFreq.index = repelem(1:size(ds.freq,1),nnz(ds.index==1));
            
            %-COMPUTE RDM
            %------------------------------------------------------------------
            %..................................................................
            params.metric = mvpa.cosmo.glm.metric.meg;
            params.center_data = mvpa.cosmo.glm.chanscaling;
            params.pseudo = mvpa.cosmo.glm.npseudo;
            params.perm = mvpa.cosmo.glm.nperm;
            params.cv = mvpa.cosmo.glm.cv;
            params.mnn = mvpa.cosmo.glm.mnn;
            %..................................................................
            mode = strcat(params.metric,'-',params.mnn,'-',params.cv,'-','meg');
            ds_dsm{f} = MeasureNeuralRDM(dsFreq,mode,params);

            Subject_RDM.(subj).(sensor){f} = ds_dsm{f};
            rdm_sum{f} = rdm_sum{f} + ds_dsm{f};
        end
        
    end
    
    %-Average RDM
    for f=1:meg.convo.frex_num
        Average_RDM.(sensor){f} = rdm_sum{f}/nsubjs;
    end
    
end

%-Save Subject RDM
%..........................................................................
for s = 1:nsubjs
    subj = meg.subject.list{s};
    
    str.file = eval(mypath.file.subj.meg.tfreq.rdm);
    str.folder = eval(mypath.folder.result.meg.rdm);
    fn.result.rdm = strcat(str.folder,'/',str.file);
    RDM = Subject_RDM.(subj);
    save(fn.result.rdm,'RDM');
end
%..........................................................................

%-Save Average RDM
%..........................................................................
str.file = eval(mypath.file.data.rdm.avg.tfreq.meg);
str.folder = eval(mypath.folder.data.rdm);
fn.result.rdm = strcat(str.folder,'/',str.file);
RDM = Average_RDM;
RDM.info.subjects = meg.subject.list;
save(fn.result.rdm,'RDM');
%..........................................................................

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-Group-Average MEG-Source RDM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(param.rdm.type,'MEG-Source')

%-IMPORT FILES AND INFORMATION
%--------------------------------------------------------------------------
%-Target and Chunk
targets = mvpa.cosmo.glm.source.targets;
chunks = mvpa.cosmo.glm.source.chunks;

nsubjs = length(meg.subject.list);
nmasks = length(mvpa.cosmo.glm.mask.name);

%-ANALYSIS
%--------------------------------------------------------------------------
%-Subject RDM
for s = 1:nsubjs
    subj = meg.subject.list{s};
    
    %-Load data
    data = load_data_meg(mypath,mvpa,mri,meg,subj);

    for m = 1:nmasks
        %-Load data set
        ds_full = read_surf_source(mypath,mri,meg,mvpa,subj,m,data);
        ds_full.sa.targets = targets;
        ds_full.sa.chunks = chunks;
        
        %-Remove constant features
        ds_full = cosmo_remove_useless_data(ds_full);
        
        %-Split data by time (chunks)
        ds_splits = cosmo_split(ds_full,'chunks');
        clear ds_full;

        for t = 1:meg.epoch.ntimes
            %-COMPUTE RDM
            %------------------------------------------------------------------
            %..................................................................
            params.metric = mvpa.cosmo.glm.metric.mri;
            params.center_data = mvpa.cosmo.glm.voxscaling;
            params.pseudo = 0;
            params.perm = 0;
            params.cv = 'none';
            params.mnn = 'none';
            %..................................................................
            mode = strcat(params.metric,'-',params.mnn,'-',params.cv,'-','mri');
            ds_dsm(:,:,t) = MeasureNeuralRDM(ds_splits{t},mode,params);
        end
        
        Subject_RDM.(subj).(mvpa.cosmo.glm.mask.name{m}) = ds_dsm;
    end
end

%-Average RDM
rdm_sum = zeros();

for m = 1:nmasks
    for s = 1:nsubjs
        subj = meg.subject.list{s};
        ds_dsm = Subject_RDM.(subj).(mvpa.cosmo.glm.mask.name{m});
        rdm_sum = rdm_sum + ds_dsm;
    end
    
    Average_RDM.(mvpa.cosmo.glm.mask.name{m}) = rdm_sum/nsubjs;
end

%-Save Subject RDM
%..........................................................................
for s = 1:nsubjs
    subj = meg.subject.list{s};
    
    str.file = eval(mypath.file.subj.meg.sources.rdm);
    str.folder = eval(mypath.folder.result.meg.rdm);
    fn.result.rdm = strcat(str.folder,'/',str.file);
    RDM = Subject_RDM.(subj);
    save(fn.result.rdm,'RDM');
end
%..........................................................................

%-Save Average RDM
%..........................................................................
str.file = eval(mypath.file.data.rdm.avg.sources.meg);
str.folder = eval(mypath.folder.data.rdm);
fn.result.rdm = strcat(str.folder,'/',str.file);
RDM = Average_RDM;
RDM.info.subjects = meg.subject.list;
save(fn.result.rdm,'RDM');
%..........................................................................

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-Feature Reweighted Group-Average MRI RDM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(param.rdm.type,'FR-MRI')

%-CHECK EXTERNALS
%--------------------------------------------------------------------------
cosmo_check_external('surfing');
cosmo_check_external('afni');
    
%-IMPORT FILES AND INFORMATION
%--------------------------------------------------------------------------
%-Target and Chunk
targets = mvpa.cosmo.glm.targets;

nsubjs = length(mri.subject.list);
nmasks = length(mvpa.cosmo.glm.mask.name);
sensors.type = {'grad'};
sensors.name = {'Gradiometer','Magnetometer'};

%-Load Group-Average MEG RDM
str.file.meg = eval(mypath.file.data.rdm.avg.meg);
str.folder = eval(mypath.folder.data.rdm);
load(strcat(str.folder,'\',str.file.meg),'RDM');

%-LOAD DATA
%==========================================================================
data = load_data_mri(mypath,mvpa,mri);

%-ANALYSIS
%==========================================================================
for sen=1:length(sensors.type)
    sensor = sensors.type{sen};
    
    for m = 1:nmasks
        msk = mvpa.cosmo.glm.mask.name{m};
        rdm_sum = zeros();

        for s = 1:nsubjs
            subj = mri.subject.list{s};

            %-Load data set
            ds_full = read_surf(mypath,mri,mvpa,subj,m,data);
            ds_full.sa.targets = targets;

            %-COMPUTE RDM
            %------------------------------------------------------------------
            %-Compute average for each unique target, so that the datasets
            %have one sample for each target
            ds = cosmo_fx(ds_full,@(x)mean(x,1),'targets',1);

            %-Remove constant features
            ds = cosmo_remove_useless_data(ds);

            %target = RDM.(sensor);
            target = mean(RDM.(sensor)(:,:,meg.epoch.nbaseline+1:meg.epoch.endstim),3);
            predictor = ds.samples;
            ds_dsm = MeasureFRRSA(mypath,target,predictor,'rdm');

            Subject_RDM.(subj).(sensor).(msk) = ds_dsm;
            rdm_sum = rdm_sum + ds_dsm;
        end

        %-Average RDM
        Average_RDM.(sensor).(msk) = rdm_sum/nsubjs;

    end

end

%-Save Subject RDM
%..........................................................................
for s = 1:nsubjs
    subj = mri.subject.list{s};
    
    str.file = eval(mypath.file.subj.frrsa.mri.rdm);
    str.folder = eval(mypath.folder.result.mri.rdm);
    fn.result.rdm = strcat(str.folder,'/',str.file);
    RDM = Subject_RDM.(subj);
    save(fn.result.rdm,'RDM');
end
%..........................................................................

%-Save Average RDM
%..........................................................................
str.file = eval(mypath.file.data.rdm.avg.frrsa.mri);
str.folder = eval(mypath.folder.data.rdm);
fn.result.rdm = strcat(str.folder,'/',str.file);
RDM = Average_RDM;
RDM.info.subjects = mri.subject.list;
save(fn.result.rdm,'RDM');
%..........................................................................

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-Feature Reweighted Group-Average MEG RDM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(param.rdm.type,'FR-MEG')

%-IMPORT FILES AND INFORMATION
%--------------------------------------------------------------------------
nsubjs = length(meg.subject.list);
nmasks = length(mvpa.cosmo.glm.mask.name);
sensors.type = {'grad','mag'};
sensors.name = {'Gradiometer','Magnetometer'};

%-Load Group-Average MRI RDM
str.file.mri = eval(mypath.file.data.rdm.avg.mri);
str.folder = eval(mypath.folder.data.rdm);
load(strcat(str.folder,'\',str.file.mri),'RDM');

%-ANALYSIS
%--------------------------------------------------------------------------
%-Iterate through sensors
for sen=1:length(sensors.type)
    sensor = sensors.type{sen};
        
    for m = 1:nmasks
        msk = mvpa.cosmo.glm.mask.name{m};
        rdm_sum = zeros();
        
        for s = 1:nsubjs
            subj = meg.subject.list{s};

            %-Load data set
            str.file = eval(mypath.file.subj.meg.preprocessed);
            str.folder = eval(mypath.folder.subj.meg);
            fn.in = strcat(str.folder,'/',str.file);
            ds = load(fn.in);

            %-Load data of Gradiometer or Magnetometer based on
            %sensor parameter.
            if strcmp(sensor,'mag')
                ds.trial = ds.trial(:,ds.channel_type(:,1)=='m',:);
            elseif strcmp(sensor,'grad')
                ds.trial = ds.trial(:,ds.channel_type(:,1)=='g',:);
            else
                warning("Please choose sensor type (mag or grad)!");
            end

            %-COMPUTE RDM
            %------------------------------------------------------------------
            %-Define targer
            target = RDM.(mvpa.cosmo.glm.mask.name{m});
            %-Define predictor
            X = ds.trial;
            y = single(ds.index);
            n_conditions = length(unique(y));
            for c = 1:n_conditions
                predictor(c,:,:) = mean(X(y==c,:,:),1);
            end
            
            for t = 1:meg.epoch.ntimes
                ds_dsm(:,:,t) = MeasureFRRSA(mypath,target,predictor(:,:,t),'rdm');
            end

            Subject_RDM.(subj).(sensor).(msk) = ds_dsm;
            rdm_sum = rdm_sum + ds_dsm;
        end
        
        %-Average RDM
        Average_RDM.(sensor).(msk) = rdm_sum/nsubjs;
        
    end
    
end

%-Save Subject RDM
%..........................................................................
for s = 1:nsubjs
    subj = meg.subject.list{s};
    
    str.file = eval();
    str.folder = eval(mypath.folder.result.meg.rdm);
    fn.result.rdm = strcat(str.folder,'/',str.file);
    RDM = Subject_RDM.(subj);
    save(fn.result.rdm,'RDM');
end
%..........................................................................

%-Save Average RDM
%..........................................................................
str.file = eval(mypath.file.data.rdm.avg.frrsa.meg);
str.folder = eval(mypath.folder.data.rdm);
fn.result.rdm = strcat(str.folder,'/',str.file);
RDM = Average_RDM;
RDM.info.subjects = meg.subject.list;
save(fn.result.rdm,'RDM');
%..........................................................................

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
                ds_full = ds;
            else
                ds_full = cosmo_stack({ds_full,ds});
            end
        end
    end
    
end

%__________________________________________________________________________
function ds_full = read_surf_source(mypath,mri,meg,mvpa,subj,m,data)

    ds_full = struct();
    
    activeVoxels = read_voxel(mypath,mri,mvpa,subj,m,data);
    
    for f = 1:length(mvpa.rdm.data.name)
        %-Select active voxels
        ds = cosmo_slice(data.(subj){f},logical(activeVoxels),2);

        if(isequal(ds_full,struct()))
            ds_full = ds;
        else
            ds_full = cosmo_stack({ds_full,ds});
        end
    end
    
end

%__________________________________________________________________________
function data = load_data_mri(mypath,mvpa,mri)

nsubjs = length(mri.subject.list);
nmasks = length(mvpa.cosmo.glm.mask.name);

%-Read Atlas
%--------------------------------------------------------------------------
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
%--------------------------------------------------------------------------
%-Read sources if it is necessary
if mri.mask.nvox~=0
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
end

%-Read Subjects' Data (1st-Level T-Map)
%--------------------------------------------------------------------------
for sbj = 1:nsubjs
    subj = mri.subject.list{sbj};
    for s = 1:mri.analysis.split.num
        for d=1:length(mvpa.cosmo.glm.file.(mri.task.name))
            c = mvpa.cosmo.glm.file.(mri.task.name){d};
            str.file = eval(mypath.file.subj.mri.firstLevel.dset);
            str.folder = eval(mypath.folder.subj.mri.firstLevel);
            fn.data = strcat(str.folder,'/',str.file);
            data.(subj){s,d} = cosmo_surface_dataset(fn.data);
        end
    end
end

end

%__________________________________________________________________________
function data = load_data_meg(mypath,mvpa,mri,meg,subj)

nmasks = length(mvpa.cosmo.glm.mask.name);

%-Read Atlas
%--------------------------------------------------------------------------
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
%--------------------------------------------------------------------------
%-Read sources if it is necessary
if mri.mask.nvox~=0
    str.file = eval(mypath.file.subj.mri.source.(mri.mask.voxel.type).dset);
    str.folder = eval(mypath.folder.result.mri.source.(mri.mask.voxel.type));
    fn.source = strcat(str.folder,'/',str.file);
    if strcmp(mri.mask.voxel.type,'firstLevel')
        data.source = cosmo_surface_dataset(fn.source);
    elseif strcmp(mri.mask.voxel.type,'secondLevel')
        s = find(strcmp(meg.subject.list,subj));
        data.source{s} = cosmo_surface_dataset(fn.source);
        data.source{s}.samples = data.source{s}.samples(1,:);
    end
end

%-Read Subjects' Data (Source MEG)
%--------------------------------------------------------------------------
for f = 1:length(mvpa.rdm.data.name)
    str.file = eval(mypath.file.subj.meg.source);
    str.folder = eval(mypath.folder.subj.source);
    fn.data = strcat(str.folder,'/',str.file);
    data.(subj){f} = cosmo_surface_dataset(fn.data);
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
function Average_Num_RDM = makeNumRDM(mvpa,Average_RDM,reg_out_vec,reg_out)

    if reg_out
        [~,~,resid] = regress(cosmo_squareform(Average_RDM,'tovector')',reg_out_vec);
        Average_RDM = cosmo_squareform(1-resid);
    end

    n_num = mvpa.cosmo.tf.cof.num;
    n_inst = mvpa.cosmo.tf.cof.inst;
    
    for i=1:n_num
        for j=1:n_num
            if i==j
                Average_Num_RDM(i,j) = ...
                    mean2(cosmo_squareform(Average_RDM((i-1)*n_inst+1:i*n_inst,(j-1)*n_inst+1:j*n_inst),'tovector'));
            else
                Average_Num_RDM(i,j) = ...
                    mean2(Average_RDM((i-1)*n_inst+1:i*n_inst,(j-1)*n_inst+1:j*n_inst));
            end
        end
    end

end

%__________________________________________________________________________
function [indiv,merged] = makeTF(mvpa,Average_Num_RDM)

    Average_Num_RDM = 1-Average_Num_RDM;
    
    %-Individual Number
    %......................................................................
    for n=1:mvpa.cosmo.tf.cof.num
        var_name = matlab.lang.makeValidName(mvpa.cosmo.tf.num.label{n});
        indiv.x.(var_name) = log(mvpa.cosmo.tf.num.value);
        indiv.y.(var_name) = normalize(Average_Num_RDM(n,:),'range');
        try
            indiv.f.(var_name) = ...
                fit(indiv.x.(var_name).',indiv.y.(var_name).','gauss1');
            indiv.fwhm.(var_name) = 2*sqrt(log(2))*indiv.f.(var_name).c1;
        end
    end

    %-Numbers Merged Together
    %......................................................................
    merged.x = [];
    merged.y = [];
    for n=1:mvpa.cosmo.tf.cof.num
        var_name = matlab.lang.makeValidName(mvpa.cosmo.tf.num.label{n});
        merged.x = [merged.x,indiv.x.(var_name)-log(mvpa.cosmo.tf.num.value(n))];
        merged.y = [merged.y,indiv.y.(var_name)];
    end
    
    try
        merged.f = fit(merged.x.',normalize(merged.y,'range').','gauss1');
        merged.fwhm = 2*sqrt(log(2))*merged.f.c1;
    end
end
