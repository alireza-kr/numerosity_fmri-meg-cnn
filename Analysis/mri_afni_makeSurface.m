%-Run the following code in the bash before executing the scrip
%--------------------------------------------------------------------------
%export MATLAB_SHELL="/bin/bash"
%export FREESURFER_HOME=/home/ali/FreeSurfer
%source $FREESURFER_HOME/SetUpFreeSurfer.sh
%export SUBJECTS_DIR=`pwd`
%conda activate py27

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-Create Surface using Surfing toolbox
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nsubjs = length(mri.subject.list);

for s=1:nsubjs
    subj = mri.subject.list{s};
    
    cmd.hdr = {'#!/bin/bash'};

    cmd.dir.anat = {'anat_dir=''%s'''};
    cmd.dir.functional = {'functional_dir=''%s'''};
    cmd.dir.freesurfer = {'freesurfer_dir=''%s'''};
    cmd.dir.afni = {'afni_dir=''%s'''};
    cmd.dir.surfing = {'surfing_dir=''%s'''};
    cmd.dir.pkg.surfing = {'pkg_surfing_dir=''%s'''};

    cmd.cmd = {...
    'mkdir ${afni_dir}'
    'mkdir ${surfing_dir}'

    '#python "${pkg_surfing_dir}"/prep_afni_surf.py -d ${freesurfer_dir} -e ${functional_dir}/meanfmri.nii -r ${surfing_dir} -l 16+64+141'
    'python "${pkg_surfing_dir}"/prep_afni_surf.py -d ${freesurfer_dir} -a ${anat_dir}/r_anat.nii -r ${surfing_dir} -l 16+64+141'

    'mv ${freesurfer_dir}/surf/SUMA/*.* ${afni_dir}'
    'rm -r ${freesurfer_dir}/surf/SUMA'};

    str.file = strcat(subj,'_surfing.sh');
    str.folder = mypath.folder.temp.root;
    fn.script = strcat(str.folder,'/',str.file);

    fileID = fopen(fn.script,'w');
    
    %-Header
    %----------------------------------------------------------------------
    fprintf(fileID,'%s\n',string(cmd.hdr));
    fprintf(fileID,'\n');
    
    %-Anatomical Path
    %----------------------------------------------------------------------
    fprintf(fileID,string(cmd.dir.anat),eval(mypath.folder.subj.mri.anatomical));
    fprintf(fileID,'\n');
    
    %-Functionl Path
    %----------------------------------------------------------------------
    fprintf(fileID,string(cmd.dir.functional),eval(mypath.folder.subj.mri.functional));
    fprintf(fileID,'\n');
    
    %-FreeSurfer Path
    %----------------------------------------------------------------------
    fprintf(fileID,string(cmd.dir.freesurfer),eval(mypath.folder.subj.surface.freesurfer));
    fprintf(fileID,'\n');
    
    %-AFNI Path
    %----------------------------------------------------------------------
    fprintf(fileID,string(cmd.dir.afni),eval(mypath.folder.subj.surface.afni));
    fprintf(fileID,'\n');
    
    %-Surfing Path
    %----------------------------------------------------------------------
    fprintf(fileID,string(cmd.dir.surfing),eval(mypath.folder.subj.surface.surfing));
    fprintf(fileID,'\n');
    
    %-Surfing Package
    %----------------------------------------------------------------------
    fprintf(fileID,string(cmd.dir.pkg.surfing),mypath.folder.code.surfing);
    fprintf(fileID,'\n');
    
    %-Main Script
    %----------------------------------------------------------------------
    fprintf(fileID,'%s\n',string(cmd.cmd));
    fclose(fileID);

    %-Run the command
    [status,cmdout] = system(sprintf('sh %s',fn.script),'-echo');

    %-Delete the script
    delete(fn.script);

end
