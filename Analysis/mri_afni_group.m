%https://andysbrainbook.readthedocs.io/en/latest/AFNI/AFNI_Short_Course/SUMA/SUMA_04_GroupAnalysisOnTheSurface.html

%__________________________________________________________________________
function mri_afni_group(mypath,mri,mvpa,mode)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-ANALYSIS MODE - GLM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(mode,'glm-surface-voxel')
    
for p=1:length(mvpa.rdm.predictor.name)

    fn.in = "";
    fn.out = "";
    
    for s=1:length(mri.subject.list)
        subj = mri.subject.list{s};
        
        str.file = eval(mypath.file.subj.mri.glm.searchlight);
        str.folder = eval(mypath.folder.result.mri.searchlight);
        fn.in = [fn.in,strcat(str.folder,'/',str.file)];
    end
    
    str.file = eval(mypath.file.result.voxel.mri.glm);
    str.folder = eval(mypath.folder.result.mri.group);
    fn.out = strcat(str.folder,'/',str.file);
    
    cmd.hdr = {'#!/bin/bash'};

    cmd.files.in = {'inputs=''%s'''};
    cmd.files.out = {'output=''%s'''};

    cmd.cmd = {...
    '3dttest++ -prefix ${output} \'
    '  -setA ${inputs}'

    'done'};

    str.file = strcat('afni_group.sh');
    str.folder = mypath.folder.temp.root;
    fn.script = strcat(str.folder,'/',str.file);

    fileID = fopen(fn.script,'w');
    fprintf(fileID,'%s\n',string(cmd.hdr));
    fprintf(fileID,'\n');
    fprintf(fileID,string(cmd.files.in),strjoin(fn.in));
    fprintf(fileID,'\n');
    fprintf(fileID,string(cmd.files.out),fn.out);
    fprintf(fileID,'\n');
    fprintf(fileID,'%s\n',string(cmd.cmd));
    fclose(fileID);
    
    [status,cmdout] = system(sprintf('sh %s',fn.script),'-echo');

    %-Delete the script
    str.file = strcat('afni_group.sh');
    str.folder = mypath.folder.temp.root;
    fn.script = strcat(str.folder,'/',str.file);
    delete(fn.script);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-ANALYSIS MODE - SECOND-LEVEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(mode,'second-level-surface-voxel')

fn.in = "";
fn.out = ""; 
    
for s=1:length(mri.subject.list)
    subj = mri.subject.list{s};

    str.file = eval(mypath.file.subj.mri.secondLevel.dset);
    str.folder = eval(mypath.folder.subj.mri.firstLevel);
    fn.in = [fn.in,strcat(str.folder,'/',str.file)];
end

str.file = eval(mypath.file.result.voxel.mri.secondLevel);
str.folder = eval(mypath.folder.result.mri.group);
fn.out = strcat(str.folder,'/',str.file);

cmd.hdr = {'#!/bin/bash'};

cmd.files.in = {'inputs=''%s'''};
cmd.files.out = {'output=''%s'''};

cmd.cmd = {...
'3dttest++ -prefix ${output} \'
'  -setA ${inputs}'

'done'};

str.file = strcat('afni_group.sh');
str.folder = mypath.folder.temp.root;
fn.script = strcat(str.folder,'/',str.file);

fileID = fopen(fn.script,'w');
fprintf(fileID,'%s\n',string(cmd.hdr));
fprintf(fileID,'\n');
fprintf(fileID,string(cmd.files.in),strjoin(fn.in));
fprintf(fileID,'\n');
fprintf(fileID,string(cmd.files.out),fn.out);
fprintf(fileID,'\n');
fprintf(fileID,'%s\n',string(cmd.cmd));
fclose(fileID);

[status,cmdout] = system(sprintf('sh %s',fn.script),'-echo');

%-Delete the script
str.file = strcat('afni_group.sh');
str.folder = mypath.folder.temp.root;
fn.script = strcat(str.folder,'/',str.file);
delete(fn.script);
    
end

end
