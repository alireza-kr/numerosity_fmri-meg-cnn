%https://brainhack-princeton.github.io/handbook/content_pages/05-03-surfaceBased.html
%https://afni.nimh.nih.gov/afni/community/board/read.php?1,79862,79862

%-Studies smoothed surface
%--------------------------------------------------------------------------
%https://www.nature.com/articles/s41598-020-62832-z
%https://elifesciences.org/articles/32962
%https://elifesciences.org/articles/66276
%https://www.pnas.org/content/118/25/e2026099118

%-Run the following code in the bash before executing the scrip
%--------------------------------------------------------------------------
%export MATLAB_SHELL="/bin/bash"
%conda activate py27

%__________________________________________________________________________
function mri_afni_makeSmooth(mypath,mri,mvpa,smooth)

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
            
            fn.data = eval(mypath.file.subj.mri.firstLevel.dset);
            files.in.con = [files.in.con,fn.data];
            files.out.con = [files.out.con,strrep(fn.data,'fwhm0',strcat('fwhm',num2str(smooth)))];
        end

        %-spmT files
        for d=1:length(mvpa.cosmo.surface.file.(mri.task.name).spmt)
            c = mvpa.cosmo.surface.file.(mri.task.name).spmt{d};
            mri.file.type = 'spmt';
            
            fn.data = eval(mypath.file.subj.mri.firstLevel.dset);
            files.in.spmt = [files.in.spmt,fn.data];
            files.out.spmt = [files.out.spmt,strrep(fn.data,'fwhm0',strcat('fwhm',num2str(smooth)))];
        end
    end
end

for sub=1:length(mri.subject.list)
    subj = mri.subject.list{sub};
    
    cmd.hdr = {'#!/bin/bash'};

    cmd.smooth = {'smooth=%g'};
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

    'SurfSmooth -spec ${surfing_dir}/mh_ico141_al.spec \'
    '  -surf_A ico141_mh.smoothwm_al.asc \'
    '  -surf_B ico141_mh.pial_al.asc \'
    '  -input ${firstlevel_dir}/${con_in_files[$i]} \'
    '  -met HEAT_07 \'
    '  -target_fwhm ${smooth} \'
    '  -Niter -1 \'
    '  -no_detrend_master \'
    '  -output ${firstlevel_dir}/${con_out_files[$i]}'

    'done'};


    cmd.cmd.spmt = {...
    'spmt_in_files=( $spmt_in_files )'
    'spmt_out_files=( $spmt_out_files )'
    
    'len=${#spmt_in_files[@]}'
    
    'for ((i=0; i<$len; i++))'
    'do'

    'SurfSmooth -spec ${surfing_dir}/mh_ico141_al.spec \'
    '  -surf_A ico141_mh.smoothwm_al.asc \'
    '  -surf_B ico141_mh.pial_al.asc \'
    '  -input ${firstlevel_dir}/${spmt_in_files[$i]} \'
    '  -met HEAT_07 \'
    '  -target_fwhm ${smooth} \'
    '  -Niter -1 \'
    '  -no_detrend_master \'
    '  -output ${firstlevel_dir}/${spmt_out_files[$i]}'

    'done'};

    str.file = strcat(subj,'_afni_smooth.sh');
    str.folder = mypath.folder.temp.root;
    fn.script = strcat(str.folder,'/',str.file);

    fileID = fopen(fn.script,'w');
    
    %-Header
    %----------------------------------------------------------------------
    fprintf(fileID,'%s\n',string(cmd.hdr));
    fprintf(fileID,'\n');
    
    %-Smoothing Level
    %----------------------------------------------------------------------
    fprintf(fileID,string(cmd.smooth),smooth);
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

    %-Run the command
    [status,cmdout] = system(sprintf('/bin/bash %s',fn.script),'-echo');

    %-Delete the script
    str.file = strcat(subj,'_afni_smooth.sh');
    str.folder = mypath.folder.temp.root;
    fn.script = strcat(str.folder,'/',str.file);
    delete(fn.script);

    %-Delete smrec files
    str.file = '*.smrec';
    str.folder = eval(mypath.folder.subj.mri.firstLevel);
    fn.smrec = strcat(str.folder,'/',str.file);
    delete(fn.smrec);
    
end
 
end
