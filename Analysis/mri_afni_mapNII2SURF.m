%__________________________________________________________________________
function mri_afni_mapNII2SURF(mypath,mri,mvpa)

for sub = 1:length(mri.subject.list)
    subj = mri.subject.list{sub};
    
    %-con files
    files.in.con = {};
    files.out.con = {};
    %-spmT files
    files.in.spmt = {};
    files.out.spmt = {};
    for s = 1:mri.analysis.split.num
        %-con files
        for d=1:length(mvpa.cosmo.surface.file.(mri.task.name).con)
            c = mvpa.cosmo.surface.file.(mri.task.name).con{d};
            mri.file.type = 'con';
            
            fn.data.in = eval(mypath.file.subj.mri.firstLevel.nii);
            fn.data.out = eval(mypath.file.subj.mri.firstLevel.dset);
            files.in.con = [files.in.con,fn.data.in];
            files.out.con = [files.out.con,fn.data.out];
        end

        %-spmT files
        for d=1:length(mvpa.cosmo.surface.file.(mri.task.name).spmt)
            c = mvpa.cosmo.surface.file.(mri.task.name).spmt{d};
            mri.file.type = 'spmt';
            
            fn.data.in = eval(mypath.file.subj.mri.firstLevel.nii);
            fn.data.out = eval(mypath.file.subj.mri.firstLevel.dset);
            files.in.spmt = [files.in.spmt,fn.data.in];
            files.out.spmt = [files.out.spmt,fn.data.out];
        end
    end

    cmd.hdr = {'#!/bin/bash'};

    cmd.dir.anat = {'anat_dir=''%s'''};
    cmd.dir.firstLevel = {'firstlevel_dir=''%s'''};
    cmd.dir.surfing = {'surfing_dir=''%s'''};
    cmd.files.in.con = {'con_in_files=''%s'''};
    cmd.files.out.con = {'con_out_files=''%s'''};
    cmd.files.in.spmt = {'spmt_in_files=''%s'''};
    cmd.files.out.spmt = {'spmt_out_files=''%s'''};

    cmd.cmd.con = {...
    'con_in_files=( $con_in_files )'
    'con_out_files=( $con_out_files )'
    
    'len=${#con_in_files[@]}'

    'for ((i=0; i<$len; i++))'
    'do'

    '3dVol2Surf -spec ${surfing_dir}/mh_ico141_al.spec \'
    '  -surf_A ico141_mh.smoothwm_al.asc \'
    '  -surf_B ico141_mh.pial_al.asc \'
    '  -sv ${anat_dir}/r_anat.nii \'
    '  -grid_parent ${firstlevel_dir}/${con_in_files[$i]} \'
    '  -f_steps 10 \'
    '  -f_index voxels \'
    '  -f_p1_fr -0.1 \'
    '  -f_pn_fr 0.2 \'
    '  -map_func ave \'
    '  -out_niml ${firstlevel_dir}/${con_out_files[$i]}'

    'done'};

    cmd.cmd.spmt = {...
    'spmt_in_files=( $spmt_in_files )'
    'spmt_out_files=( $spmt_out_files )'
    
    'len=${#spmt_in_files[@]}'
    
    'for ((i=0; i<$len; i++))'
    'do'

    '3dVol2Surf -spec ${surfing_dir}/mh_ico141_al.spec \'
    '  -surf_A ico141_mh.smoothwm_al.asc \'
    '  -surf_B ico141_mh.pial_al.asc \'
    '  -sv ${anat_dir}/r_anat.nii \'
    '  -grid_parent ${firstlevel_dir}/${spmt_in_files[$i]} \'
    '  -f_steps 10 \'
    '  -f_index voxels \'
    '  -f_p1_fr -0.1 \'
    '  -f_pn_fr 0.2 \'
    '  -map_func ave \'
    '  -out_niml ${firstlevel_dir}/${spmt_out_files[$i]}'

    'done'};

    str.file = strcat(subj,'_nii2surf.sh');
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
    
    %-1st-Level Path
    %----------------------------------------------------------------------
    fprintf(fileID,string(cmd.dir.firstLevel),eval(mypath.folder.subj.mri.firstLevel));
    fprintf(fileID,'\n');
    
    %-Surfing Path
    %----------------------------------------------------------------------
    fprintf(fileID,string(cmd.dir.surfing),eval(mypath.folder.subj.surface.surfing));
    fprintf(fileID,'\n');

    %-Input Files (con)
    %----------------------------------------------------------------------
    fprintf(fileID,string(cmd.files.in.con),strjoin(files.in.con));
    fprintf(fileID,'\n');
    
    %-Output Files (con)
    %----------------------------------------------------------------------
    fprintf(fileID,string(cmd.files.out.con),strjoin(files.out.con));
    fprintf(fileID,'\n');
    
    %-Input Files (spmT)
    %----------------------------------------------------------------------
    fprintf(fileID,string(cmd.files.in.spmt),strjoin(files.in.spmt));
    fprintf(fileID,'\n');
    
    %-Output Files (spmT)
    %----------------------------------------------------------------------
    fprintf(fileID,string(cmd.files.out.spmt),strjoin(files.out.spmt));
    fprintf(fileID,'\n');
    
    %-Main Script (for con files)
    %----------------------------------------------------------------------
    fprintf(fileID,'%s\n',string(cmd.cmd.con));
    
    %-Main Script (for spmT files)
    %----------------------------------------------------------------------
    fprintf(fileID,'%s\n',string(cmd.cmd.spmt));
    
    fclose(fileID);
    
    [status,cmdout] = system(sprintf('/bin/bash %s',fn.script),'-echo');

    %-Delete the script
    str.file = strcat(subj,'_nii2surf.sh');
    str.folder = mypath.folder.temp.root;
    fn.script = strcat(str.folder,'/',str.file);
    delete(fn.script);

end

end
