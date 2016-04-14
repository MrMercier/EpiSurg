% load MRI header
MRIhdr=MRIread(fullfile(getenv('SUBJECTS_DIR'),subj,'/mri/orig_1.mgh'),true);

% the transform matrices
VOX2RAS=[-1 0 0 128; 0 0 -1 128; 0 -1 0 128; 0 0 0 1];
TalXFM=freesurfer_read_talxfm(fullfile(getenv('SUBJECTS_DIR'),subj,'/mri/transforms/talairach.xfm'));
Norig=MRIhdr.vox2ras;
Torig=MRIhdr.tkrvox2ras;

% put channel coords in RAS space
chanpos_vox=the_coord(the_mask_depth,:);
chanpos_ras=(VOX2RAS*[chanpos_vox'; ones(1, size(chanpos_vox,1))])';
chanpos_ras(:,4)=[];

% the actual transform!
% code taken from http://surfer.nmr.mgh.harvard.edu/fswiki/CoordinateSystems
% Transforms within a subject's anatomical space
% 2. I have an RAS point on the surface (tkrR tkrA tkrS) ("Vertex RAS" from tksurfer) and want to compute the MNI305 RAS that corresponds to this point:
% MNI305RAS = TalXFM*Norig*inv(Torig)*[tkrR tkrA tkrS 1]'
% TalXFM: subject/orig/transforms/talairach.xfm Norig: mri_info --vox2ras orig.mgz Torig: mri_info --vox2ras-tkr orig.mgz
MNI305RAS_chan=(TalXFM*Norig*inv(Torig)*[chanpos_ras'; ones(1, size(chanpos_ras,1))])';
MNI305RAS_chan(:,4)=[];

% and now back into PIALVOX space
MNI305VOX_chan=(VOX2RAS\[MNI305RAS_chan'; ones(1, size(MNI305RAS_chan,1))])';
MNI305VOX_chan(:,4)=[];


