%-Run the following package
%--------------------------------------------------------------------------
%conda create -n stc2gii python=3.8
%conda activate stc2gii
%pip install stc2gii_hack

%-Run the following code in the bash before executing the scrip
%--------------------------------------------------------------------------
%conda activate stc2gii

%https://mne.discourse.group/t/best-practice-to-get-niml-gifti-or-nifti-from-stc/5767/3
%https://github.com/jbteves/stc2gii_hack

%__________________________________________________________________________
function meg_afni_mapSTC2SURF(mypath,mri,meg,mvpa,smooth)

for sub = 1:length(meg.subject.list)
    subj = meg.subject.list{sub};
    
    %-Copy the std.141 files to the temporary folder
    %----------------------------------------------------------------------
    try
        copyfile...
            (strcat(eval(mypath.folder.subj.surface.freesurfer),'/surf/SUMA/','std.141.lh.white.gii'),...
            eval(mypath.folder.temp.root));
        copyfile...
            (strcat(eval(mypath.folder.subj.surface.freesurfer),'/surf/SUMA/','std.141.rh.white.gii'),...
            eval(mypath.folder.temp.root));
        copyfile...
            (strcat(eval(mypath.folder.subj.surface.freesurfer),'/surf/SUMA/','std.141.lh.pial.gii'),...
            eval(mypath.folder.temp.root));
        copyfile...
            (strcat(eval(mypath.folder.subj.surface.freesurfer),'/surf/SUMA/','std.141.rh.pial.gii'),...
            eval(mypath.folder.temp.root));
        copyfile...
            (strcat(eval(mypath.folder.subj.surface.freesurfer),'/surf/SUMA/','std.141.lh.smoothwm.gii'),...
            eval(mypath.folder.temp.root));
        copyfile...
            (strcat(eval(mypath.folder.subj.surface.freesurfer),'/surf/SUMA/','std.141.rh.smoothwm.gii'),...
            eval(mypath.folder.temp.root));
        copyfile...
            (strcat(eval(mypath.folder.subj.surface.freesurfer),'/surf/SUMA/','std.141.freesurfer_lh.spec'),...
            eval(mypath.folder.temp.root));
        copyfile...
            (strcat(eval(mypath.folder.subj.surface.freesurfer),'/surf/SUMA/','std.141.freesurfer_rh.spec'),...
            eval(mypath.folder.temp.root));
    catch

        %..............................................................
        cmd.hdr = {'#!/bin/bash'};

        cmd.dir.surf = {'folder_surf=''%s'''};

        cmd.cmd = {...
        '@SUMA_Make_Spec_FS -fspath $folder_surf -sid freesurfer -GIFTI'};

        str.file = strcat(subj,'_surf.sh');
        str.folder = eval(mypath.folder.temp.root);
        fn.script = strcat(str.folder,'/',str.file);

        fileID = fopen(fn.script,'w');

        %-Header
        %--------------------------------------------------------------
        fprintf(fileID,'%s\n',string(cmd.hdr));
        fprintf(fileID,'\n');
        %..............................................................

        %-Surf Path (inside FreeSurfer Folder)
        %--------------------------------------------------------------
        fprintf(fileID,string(cmd.dir.surf),strcat(eval(mypath.folder.subj.surface.freesurfer),'/surf/'));
        fprintf(fileID,'\n');

        %-Main Script
        %----------------------------------------------------------------------
        fprintf(fileID,'%s\n',string(cmd.cmd));
        fprintf(fileID,'\n');

        fclose(fileID);
        
        %-Run the command
        %https://www.mathworks.com/matlabcentral/answers/157199-run-a-bash-script-from-matlab
        [status,cmdout] = system(sprintf('sh %s',fn.script),'-echo');

        %-Delete the script
        delete(fn.script);

        copyfile...
            (strcat(eval(mypath.folder.subj.surface.freesurfer),'/surf/SUMA/','std.141.lh.white.gii'),...
            eval(mypath.folder.temp.root));
        copyfile...
            (strcat(eval(mypath.folder.subj.surface.freesurfer),'/surf/SUMA/','std.141.rh.white.gii'),...
            eval(mypath.folder.temp.root));
        copyfile...
            (strcat(eval(mypath.folder.subj.surface.freesurfer),'/surf/SUMA/','std.141.lh.pial.gii'),...
            eval(mypath.folder.temp.root));
        copyfile...
            (strcat(eval(mypath.folder.subj.surface.freesurfer),'/surf/SUMA/','std.141.rh.pial.gii'),...
            eval(mypath.folder.temp.root));
        copyfile...
            (strcat(eval(mypath.folder.subj.surface.freesurfer),'/surf/SUMA/','std.141.lh.smoothwm.gii'),...
            eval(mypath.folder.temp.root));
        copyfile...
            (strcat(eval(mypath.folder.subj.surface.freesurfer),'/surf/SUMA/','std.141.rh.smoothwm.gii'),...
            eval(mypath.folder.temp.root));
        copyfile...
            (strcat(eval(mypath.folder.subj.surface.freesurfer),'/surf/SUMA/','std.141.freesurfer_lh.spec'),...
            eval(mypath.folder.temp.root));
        copyfile...
            (strcat(eval(mypath.folder.subj.surface.freesurfer),'/surf/SUMA/','std.141.freesurfer_rh.spec'),...
            eval(mypath.folder.temp.root));
    end
    
    %-Copy the STC files to the temporary folder
    %----------------------------------------------------------------------
    for f = 1:length(mvpa.rdm.data.name)
        % Left STC
        copyfile...
            (strcat(eval(mypath.folder.subj.source),'/',eval(mypath.file.subj.meg.stc.lh)),...
            eval(mypath.folder.temp.root));
        % Right STC
        copyfile...
            (strcat(eval(mypath.folder.subj.source),'/',eval(mypath.file.subj.meg.stc.rh)),...
            eval(mypath.folder.temp.root));
    end
    
    %-Do the conversion (STC to GII)
    %----------------------------------------------------------------------
    clear str;
    clear cmd;
    
    files.in.stc.lh = {};
    files.in.stc.rh = {};
    files.in.src = {};
    files.out = {};
    
    for f = 1:length(mvpa.rdm.data.name)
        files.in.stc.lh = [files.in.stc.lh,eval(mypath.file.subj.meg.stc.lh)];
        files.in.stc.rh = [files.in.stc.rh,eval(mypath.file.subj.meg.stc.rh)];
        files.in.src = [files.in.src,eval(mypath.file.subj.meg.src)];
        files.out = [files.out,strcat(mvpa.rdm.data.name{f})];
    end
    
    cmd.hdr = {'#!/bin/bash'};
    cmd.dir.source = {'source_dir=''%s'''};
    cmd.dir.temp = {'temp_dir=''%s'''};
    cmd.files.in.stc.lh = {'lh_in_files=''%s'''};
    cmd.files.in.stc.rh = {'rh_in_files=''%s'''};
    cmd.files.in.src = {'src_in_files=''%s'''};
    cmd.files.out = {'out_files=''%s'''};
    
    cmd.cmd = {...
    'lh_in_files=( $lh_in_files )'
    'rh_in_files=( $rh_in_files )'
    'src_in_files=( $src_in_files )'
    'out_files=( $out_files )'
    
    'len=${#out_files[@]}'

    'for ((i=0; i<$len; i++))'
    'do'
    
    'stc2gii_hack \'
    '    ${source_dir}/${src_in_files[$i]} \'
    '    ${temp_dir}/${lh_in_files[$i]} \'
    '    ${temp_dir}/${rh_in_files[$i]} \'
    '    ${temp_dir}/${out_files[$i]}'
    
    'done'};
    
    str.file = strcat(subj,'_stc2gii.sh');
    str.folder = eval(mypath.folder.temp.root);
    fn.script = strcat(str.folder,'/',str.file);
    
    fileID = fopen(fn.script,'w');
    
    %-Header
    %----------------------------------------------------------------------
    fprintf(fileID,'%s\n',string(cmd.hdr));
    fprintf(fileID,'\n');
    
    %-Source Path
    %----------------------------------------------------------------------
    fprintf(fileID,string(cmd.dir.source),eval(mypath.folder.subj.source));
    fprintf(fileID,'\n');
    
    %-Temporary Path
    %----------------------------------------------------------------------
    fprintf(fileID,string(cmd.dir.temp),eval(mypath.folder.temp.root));
    fprintf(fileID,'\n');
    
    %-Input Files (Left)
    %----------------------------------------------------------------------
    fprintf(fileID,string(cmd.files.in.stc.lh),strjoin(files.in.stc.lh));
    fprintf(fileID,'\n');
    
    %-Input Files (Right)
    %----------------------------------------------------------------------
    fprintf(fileID,string(cmd.files.in.stc.rh),strjoin(files.in.stc.rh));
    fprintf(fileID,'\n');
    
    %-Input Files (SRC)
    %----------------------------------------------------------------------
    fprintf(fileID,string(cmd.files.in.src),strjoin(files.in.src));
    fprintf(fileID,'\n');
    
    %-Output Files
    %----------------------------------------------------------------------
    fprintf(fileID,string(cmd.files.out),strjoin(files.out));
    fprintf(fileID,'\n');
    
    %-Main Script
    %----------------------------------------------------------------------
    fprintf(fileID,'%s\n',string(cmd.cmd));
    
    fclose(fileID);
    
    [status,cmdout] = system(sprintf('/bin/bash %s',fn.script),'-echo');

    %-Delete the script
    str.file = strcat(subj,'_stc2gii.sh');
    str.folder = eval(mypath.folder.temp.root);
    fn.script = strcat(str.folder,'/',str.file);
    delete(fn.script);
    
    %-Do the conversion (GII to NIML.DSET)
    %----------------------------------------------------------------------
    clear str;
    clear files;
    clear cmd;
    
    %-GII files (lh)
    files.in.lh = {};
    files.time.lh = {};
    files.centered.lh = {};
    files.out.lh = {};
    %-GII files (rh)
    files.in.rh = {};
    files.time.rh = {};
    files.centered.rh = {};
    files.out.rh = {};
    
    for f = 1:length(mvpa.rdm.data.name)
        files.in.lh = [files.in.lh,strcat(mvpa.rdm.data.name{f},'-lh.gii')];
        files.time.lh = [files.time.lh,strcat(mvpa.rdm.data.name{f},'-lh.time.gii')];
        files.centered.lh = [files.centered.lh,strcat(mvpa.rdm.data.name{f},'-lh_centered.gii')];
        files.out.lh = [files.out.lh,strcat(mvpa.rdm.data.name{f},'-lh.niml.dset')];
        files.in.rh = [files.in.rh,strcat(mvpa.rdm.data.name{f},'-rh.gii')];
        files.time.rh = [files.time.rh,strcat(mvpa.rdm.data.name{f},'-rh.time.gii')];
        files.centered.rh = [files.centered.rh,strcat(mvpa.rdm.data.name{f},'-rh_centered.gii')];
        files.out.rh = [files.out.rh,strcat(mvpa.rdm.data.name{f},'-rh.niml.dset')];
    end
    
    cmd.hdr = {'#!/bin/bash'};
    cmd.dir.temp = {'temp_dir=''%s'''};
    cmd.files.in.lh = {'lh_in_files=''%s'''};
    cmd.files.time.lh = {'lh_time_files=''%s'''};
    cmd.files.centered.lh = {'lh_centered_files=''%s'''};
    cmd.files.out.lh = {'lh_out_files=''%s'''};
    cmd.files.in.rh = {'rh_in_files=''%s'''};
    cmd.files.time.rh = {'rh_time_files=''%s'''};
    cmd.files.centered.rh = {'rh_centered_files=''%s'''};
    cmd.files.out.rh = {'rh_out_files=''%s'''};
    
    cmd.cmd.lh = {...
    'cd ${temp_dir}'
    'lh_in_files=( $lh_in_files )'
    'lh_time_files=( $lh_time_files )'
    'lh_centered_files=( $lh_centered_files )'
    'lh_out_files=( $lh_out_files )'
    
    'len=${#lh_in_files[@]}'

    'for ((i=0; i<$len; i++))'
    'do'
    
    'SurfaceMetrics -i ${temp_dir}/${lh_in_files[$i]} -coords'

    'centerX=`3dBrickStat -mean ${temp_dir}/${lh_in_files[$i]}.coord.1D.dset[1]`'
    'centerY=`3dBrickStat -mean ${temp_dir}/${lh_in_files[$i]}.coord.1D.dset[2]`'
    'centerZ=`3dBrickStat -mean ${temp_dir}/${lh_in_files[$i]}.coord.1D.dset[3]`'
    'echo "Center for your dataset: ${centerX} ${centerY} ${centerZ}"'

    'SurfaceMetrics -i ${temp_dir}/std.141.lh.white.gii -coords'
    'stdCenterX=`3dBrickStat -mean ${temp_dir}/std.141.lh.white.gii.coord.1D.dset[1]`'
    'stdCenterY=`3dBrickStat -mean ${temp_dir}/std.141.lh.white.gii.coord.1D.dset[2]`'
    'stdCenterZ=`3dBrickStat -mean ${temp_dir}/std.141.lh.white.gii.coord.1D.dset[3]`'
    'echo "Center for your std.141: ${stdCenterX} ${stdCenterY} ${stdCenterZ}"'

    'scaleX=`ccalc ${stdCenterX} - ${centerX}`'
    'scaleY=`ccalc ${stdCenterY} - ${centerY}`'
    'scaleZ=`ccalc ${stdCenterZ} - ${centerZ}`'
    'echo "Shifting..."'
    'echo "${scaleX} ${scaleY} ${scaleZ}"'
    'if [ -e center_al.1D ]; then'
    '    rm center_al.1D'
    'fi'
    'echo "1 0 0 ${scaleX}" >> center_al.1D'
    'echo "0 1 0 ${scaleY}" >> center_al.1D'
    'echo "0 0 1 ${scaleZ}" >> center_al.1D'

    'ConvertSurface -xmat_1D center_al.1D -i ${temp_dir}/${lh_in_files[$i]} -o ${temp_dir}/${lh_centered_files[$i]}'

    'SurfToSurf -i_gii ${temp_dir}/std.141.lh.white.gii \'
    '-i_gii ${temp_dir}/${lh_centered_files[$i]} \'
    '-dset ${temp_dir}/${lh_time_files[$i]} \'
    '-prefix std.141.${lh_out_files[$i]}'
    
    'done'};

    cmd.cmd.rh = {...
    'cd ${temp_dir}'
    'rh_in_files=( $rh_in_files )'
    'rh_time_files=( $rh_time_files )'
    'rh_centered_files=( $rh_centered_files )'
    'rh_out_files=( $rh_out_files )'
    
    'len=${#rh_in_files[@]}'
    
    'for ((i=0; i<$len; i++))'
    'do'
    
    'SurfaceMetrics -i ${temp_dir}/${rh_in_files[$i]} -coords'

    'centerX=`3dBrickStat -mean ${temp_dir}/${rh_in_files[$i]}.coord.1D.dset[1]`'
    'centerY=`3dBrickStat -mean ${temp_dir}/${rh_in_files[$i]}.coord.1D.dset[2]`'
    'centerZ=`3dBrickStat -mean ${temp_dir}/${rh_in_files[$i]}.coord.1D.dset[3]`'
    'echo "Center for your dataset: ${centerX} ${centerY} ${centerZ}"'

    'SurfaceMetrics -i ${temp_dir}/std.141.rh.white.gii -coords'
    'stdCenterX=`3dBrickStat -mean ${temp_dir}/std.141.rh.white.gii.coord.1D.dset[1]`'
    'stdCenterY=`3dBrickStat -mean ${temp_dir}/std.141.rh.white.gii.coord.1D.dset[2]`'
    'stdCenterZ=`3dBrickStat -mean ${temp_dir}/std.141.rh.white.gii.coord.1D.dset[3]`'
    'echo "Center for your std.141: ${stdCenterX} ${stdCenterY} ${stdCenterZ}"'

    'scaleX=`ccalc ${stdCenterX} - ${centerX}`'
    'scaleY=`ccalc ${stdCenterY} - ${centerY}`'
    'scaleZ=`ccalc ${stdCenterZ} - ${centerZ}`'
    'echo "Shifting..."'
    'echo "${scaleX} ${scaleY} ${scaleZ}"'
    'if [ -e center_al.1D ]; then'
    '    rm center_al.1D'
    'fi'
    'echo "1 0 0 ${scaleX}" >> center_al.1D'
    'echo "0 1 0 ${scaleY}" >> center_al.1D'
    'echo "0 0 1 ${scaleZ}" >> center_al.1D'

    'ConvertSurface -xmat_1D center_al.1D -i ${temp_dir}/${rh_in_files[$i]} -o ${temp_dir}/${rh_centered_files[$i]}'

    'SurfToSurf -i_gii ${temp_dir}/std.141.rh.white.gii \'
    '-i_gii ${temp_dir}/${rh_centered_files[$i]} \'
    '-dset ${temp_dir}/${rh_time_files[$i]} \'
    '-prefix std.141.${rh_out_files[$i]}'
    
    'done'};

    str.file = strcat(subj,'_gii2dset.sh');
    str.folder = eval(mypath.folder.temp.root);
    fn.script = strcat(str.folder,'/',str.file);
    
    fileID = fopen(fn.script,'w');
    
    %-Header
    %----------------------------------------------------------------------
    fprintf(fileID,'%s\n',string(cmd.hdr));
    fprintf(fileID,'\n');
    
    %-Temporary Path
    %----------------------------------------------------------------------
    fprintf(fileID,string(cmd.dir.temp),eval(mypath.folder.temp.root));
    fprintf(fileID,'\n');
    
    %-Input Files (Left)
    %----------------------------------------------------------------------
    fprintf(fileID,string(cmd.files.in.lh),strjoin(files.in.lh));
    fprintf(fileID,'\n');
    
    %-Input Files (Right)
    %----------------------------------------------------------------------
    fprintf(fileID,string(cmd.files.in.rh),strjoin(files.in.rh));
    fprintf(fileID,'\n');
    
    %-Time Files (Left)
    %----------------------------------------------------------------------
    fprintf(fileID,string(cmd.files.time.lh),strjoin(files.time.lh));
    fprintf(fileID,'\n');
    
    %-Time Files (Right)
    %----------------------------------------------------------------------
    fprintf(fileID,string(cmd.files.time.rh),strjoin(files.time.rh));
    fprintf(fileID,'\n');
    
    %-Centered Files (Left)
    %----------------------------------------------------------------------
    fprintf(fileID,string(cmd.files.centered.lh),strjoin(files.centered.lh));
    fprintf(fileID,'\n');
    
    %-Centered Files (Right)
    %----------------------------------------------------------------------
    fprintf(fileID,string(cmd.files.centered.rh),strjoin(files.centered.rh));
    fprintf(fileID,'\n');
    
    %-Output Files (Left)
    %----------------------------------------------------------------------
    fprintf(fileID,string(cmd.files.out.lh),strjoin(files.out.lh));
    fprintf(fileID,'\n');
    
    %-Output Files (Right)
    %----------------------------------------------------------------------
    fprintf(fileID,string(cmd.files.out.rh),strjoin(files.out.rh));
    fprintf(fileID,'\n');
    
    %-Main Script (Left)
    %----------------------------------------------------------------------
    fprintf(fileID,'%s\n',string(cmd.cmd.lh));
    
    %-Main Script (Right)
    %----------------------------------------------------------------------
    fprintf(fileID,'%s\n',string(cmd.cmd.rh));
    
    fclose(fileID);
    
    [status,cmdout] = system(sprintf('/bin/bash %s',fn.script),'-echo');

    %-Delete the script
    str.file = strcat(subj,'_gii2dset.sh');
    str.folder = eval(mypath.folder.temp.root);
    fn.script = strcat(str.folder,'/',str.file);
    delete(fn.script);

    %-Smooth the surfaces
    %----------------------------------------------------------------------
    clear str;
    clear files;
    clear cmd;
    
    %-DSET files (lh)
    files.in.lh = {};
    files.out.lh = {};
    %-DSET files (rh)
    files.in.rh = {};
    files.out.rh = {};
    
    for f = 1:length(mvpa.rdm.data.name)
        files.in.lh = [files.in.lh,strcat('std.141.',mvpa.rdm.data.name{f},'-lh.',mvpa.rdm.data.name{f},'-lh.time.niml.dset')];
        files.in.rh = [files.in.rh,strcat('std.141.',mvpa.rdm.data.name{f},'-rh.',mvpa.rdm.data.name{f},'-rh.time.niml.dset')];
        files.out.lh = [files.out.lh,strcat('std.141.',num2str(smooth),'.',mvpa.rdm.data.name{f},'-lh.',mvpa.rdm.data.name{f},'-lh.time.niml.dset')];
        files.out.rh = [files.out.rh,strcat('std.141.',num2str(smooth),'.',mvpa.rdm.data.name{f},'-rh.',mvpa.rdm.data.name{f},'-rh.time.niml.dset')];
    end
    
    cmd.hdr = {'#!/bin/bash'};

    cmd.smooth = {'smooth=%g'};
    cmd.dir.temp = {'temp_dir=''%s'''};
    cmd.files.in.lh = {'lh_in_files=''%s'''};
    cmd.files.in.rh = {'rh_in_files=''%s'''};
    cmd.files.out.lh = {'lh_out_files=''%s'''};
    cmd.files.out.rh = {'rh_out_files=''%s'''};

    cmd.cmd.lh = {...
    'lh_in_files=( $lh_in_files )'
    'lh_out_files=( $lh_out_files )'
    
    'len=${#lh_in_files[@]}'
    
    'for ((i=0; i<$len; i++))'
    'do'

    'SurfSmooth -spec ${temp_dir}/std.141.freesurfer_lh.spec \'
    '  -surf_A std.141.lh.smoothwm.gii \'
    '  -surf_B std.141.lh.pial.gii \'
    '  -input ${temp_dir}/${lh_in_files[$i]} \'
    '  -met HEAT_07 \'
    '  -target_fwhm ${smooth} \'
    '  -Niter -1 \'
    '  -no_detrend_master \'
    '  -output ${temp_dir}/${lh_out_files[$i]}'

    'done'};

    cmd.cmd.rh = {...
    'rh_in_files=( $rh_in_files )'
    'rh_out_files=( $rh_out_files )'
    
    'len=${#rh_in_files[@]}'
    
    'for ((i=0; i<$len; i++))'
    'do'

    'SurfSmooth -spec ${temp_dir}/std.141.freesurfer_rh.spec \'
    '  -surf_A std.141.rh.smoothwm.gii \'
    '  -surf_B std.141.rh.pial.gii \'
    '  -input ${temp_dir}/${rh_in_files[$i]} \'
    '  -met HEAT_07 \'
    '  -target_fwhm ${smooth} \'
    '  -Niter -1 \'
    '  -no_detrend_master \'
    '  -output ${temp_dir}/${rh_out_files[$i]}'

    'done'};

    str.file = strcat(subj,'_afni_smooth.sh');
    str.folder = eval(mypath.folder.temp.root);
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
    
    %-Temporary Path
    %----------------------------------------------------------------------
    fprintf(fileID,string(cmd.dir.temp),eval(mypath.folder.temp.root));
    fprintf(fileID,'\n');
    
    %-Input Files (Left)
    %----------------------------------------------------------------------
    fprintf(fileID,string(cmd.files.in.lh),strjoin(files.in.lh));
    fprintf(fileID,'\n');
    
    %-Input Files (Right)
    %----------------------------------------------------------------------
    fprintf(fileID,string(cmd.files.in.rh),strjoin(files.in.rh));
    fprintf(fileID,'\n');
    
    %-Output Files (Left)
    %----------------------------------------------------------------------
    fprintf(fileID,string(cmd.files.out.lh),strjoin(files.out.lh));
    fprintf(fileID,'\n');
    
    %-Output Files (Right)
    %----------------------------------------------------------------------
    fprintf(fileID,string(cmd.files.out.rh),strjoin(files.out.rh));
    fprintf(fileID,'\n');
    
    %-Main Script (Left)
    %----------------------------------------------------------------------
    fprintf(fileID,'%s\n',string(cmd.cmd.lh));
    
    %-Main Script (Right)
    %----------------------------------------------------------------------
    fprintf(fileID,'%s\n',string(cmd.cmd.rh));
    
    fclose(fileID);
    
    [status,cmdout] = system(sprintf('/bin/bash %s',fn.script),'-echo');

    %-Delete the script
    str.file = strcat(subj,'_afni_smooth.sh');
    str.folder = eval(mypath.folder.temp.root);
    fn.script = strcat(str.folder,'/',str.file);
    delete(fn.script);
    
    %-Concatenate left and right DSET hemisphere
    %----------------------------------------------------------------------
    %https://afni.nimh.nih.gov/afni/community/board/read.php?1,149391,149391
    %http://cosmomvpa.org/faq.html#merge-surface-data-from-two-hemispheres
    clear str;
    clear files;
    clear cmd;
    
    %-No Smoothing
    %----------------------------------------------------------------------
    mri.analysis.smoothing = 0;
    
    for f = 1:length(mvpa.rdm.data.name)
        str.file.left = strcat('std.141.',mvpa.rdm.data.name{f},'-lh.',mvpa.rdm.data.name{f},'-lh.time.niml.dset');
        str.file.right = strcat('std.141.',mvpa.rdm.data.name{f},'-rh.',mvpa.rdm.data.name{f},'-rh.time.niml.dset');
        str.folder = eval(mypath.folder.temp.root);
        fn.in.left = strcat(str.folder,'/',str.file.left);
        fn.in.right = strcat(str.folder,'/',str.file.right);
        ds_left=cosmo_surface_dataset(fn.in.left);
        ds_right=cosmo_surface_dataset(fn.in.right);

        [unused, index]=cosmo_dim_find(ds_left, 'node_indices');
        nverts_left=max(ds_left.a.fdim.values{index});

        % get the offset to set the feature attribute index later
        offset_left=numel(ds_left.a.fdim.values{index});

        % update node indices to support indexing data from two hemispheres
        node_indices=[ds_left.a.fdim.values{index}, ...
                        nverts_left+ds_right.a.fdim.values{index}];
        ds_left.a.fdim.values{index}=node_indices;
        ds_right.a.fdim.values{index}=node_indices;

        % update node indices for right hemisphere
        assert(all(ds_left.fa.node_indices<=offset_left)); % safety check
        ds_right.fa.node_indices=ds_right.fa.node_indices+offset_left;

        % merge hemisphes
        ds_left_right=cosmo_stack({ds_left,ds_right},2);

        % save the result
        str.file.out = eval(mypath.file.subj.meg.source);
        str.folder = eval(mypath.folder.subj.source);
        fn.out.dset = strcat(str.folder,'/',str.file.out);
        cosmo_map2surface(ds_left_right,fn.out.dset);
    end
    
    %-Smoothed
    %----------------------------------------------------------------------
    mri.analysis.smoothing = smooth;
    
    for f = 1:length(mvpa.rdm.data.name)
        str.file.left = strcat('std.141.',num2str(smooth),'.',mvpa.rdm.data.name{f},'-lh.',mvpa.rdm.data.name{f},'-lh.time.niml.dset');
        str.file.right = strcat('std.141.',num2str(smooth),'.',mvpa.rdm.data.name{f},'-rh.',mvpa.rdm.data.name{f},'-rh.time.niml.dset');
        str.folder = eval(mypath.folder.temp.root);
        fn.in.left = strcat(str.folder,'/',str.file.left);
        fn.in.right = strcat(str.folder,'/',str.file.right);
        ds_left=cosmo_surface_dataset(fn.in.left);
        ds_right=cosmo_surface_dataset(fn.in.right);

        [unused, index]=cosmo_dim_find(ds_left, 'node_indices');
        nverts_left=max(ds_left.a.fdim.values{index});

        % get the offset to set the feature attribute index later
        offset_left=numel(ds_left.a.fdim.values{index});

        % update node indices to support indexing data from two hemispheres
        node_indices=[ds_left.a.fdim.values{index}, ...
                        nverts_left+ds_right.a.fdim.values{index}];
        ds_left.a.fdim.values{index}=node_indices;
        ds_right.a.fdim.values{index}=node_indices;

        % update node indices for right hemisphere
        assert(all(ds_left.fa.node_indices<=offset_left)); % safety check
        ds_right.fa.node_indices=ds_right.fa.node_indices+offset_left;

        % merge hemisphes
        ds_left_right=cosmo_stack({ds_left,ds_right},2);

        % save the result
        str.file.out = eval(mypath.file.subj.meg.source);
        str.folder = eval(mypath.folder.subj.source);
        fn.out.dset = strcat(str.folder,'/',str.file.out);
        cosmo_map2surface(ds_left_right,fn.out.dset);
    end
    
    %-Delete all files in the temporary folder
    %----------------------------------------------------------------------
    delete(strcat(eval(mypath.folder.temp.root),'/*'));
    
end

end
