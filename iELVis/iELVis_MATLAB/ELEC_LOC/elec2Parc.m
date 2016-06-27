%  elec2Parc() - Returns the name of the nearest brain area according
%                to FreeSurfer pial surface-compatible parcellations.
%
% Usage:
%  >>elecParc=elec2Parc(subj,atlas);
%
% Required Input:
%  subj  -Name of the subject's freesurfer directory (full path not
%             needed)
%
% Optional Input:
%  atlas -{'DK','D','Y7','Y17'} Anatomical atlas to use:
%           'DK'=Desikan-Killiany {default}
%           'D' =Destrieux
%           'Y7'=Yeo 7-network resting state fMRI atlas
%           'Y17'=Yeo 17-network resting state fMRI atlas
%
% Output:
%   elecParc - 2D cell array containing electrode names and their
%                 assigned cortical area
%
% Example:
% >> elecParc=elec2Parc('PT001','D');
%
% Note, depth electrode anatomical locations are taken from the subject's
% FreeSurfer parcellation in aparc+aseg.mgz (regardless of the surface
% atlas specified).
%
% Author: David M. Groppe
% Jan, 2016
% Honeylab
% University of Toronto

%
% Future Work:
% -I should have a function called neat_labels.m that will format the DK
% atlas names to be more kind to the eye. It might help to incoporate that.


function elecParc=elec2Parc(subj,atlas)

if nargin<2,
    atlas='DK';
end

fsDir=getFsurfSubDir();


% Folder with surface files
surfaceFolder=fullfile(fsDir,subj,'surf');

% Folder with cortical parcellation files
labelFolder=fullfile(fsDir,subj,'label');


%% Import electrode locations
% Pial coordinates
pialFname=fullfile(fsDir,subj,'elec_recon',sprintf('%s.PIAL',subj));
pialCoordStr=csv2Cell(pialFname,' ',2);
nElec=size(pialCoordStr,1);
pialCoord=zeros(nElec,3);
for a=1:nElec,
    for b=1:3,
        pialCoord(a,b)=str2double(pialCoordStr{a,b});
    end
end

% Need to get brainmask dimensions for flipping 3rd pvox coordinate
mriFname=fullfile(fsDir,subj,'mri','brainmask.mgz');
if ~exist(mriFname,'file')
   error('File %s not found.',mriFname); 
end
mri=MRIread(mriFname);
%mri.vol is ILA (i.e., S->I, R->L, P->A)
sVol=size(mri.vol);
clear mri mriFname

% Voxel coordinates
pvoxFname=fullfile(fsDir,subj,'elec_recon',sprintf('%s.PIALVOX',subj));
pvoxCoordStr=csv2Cell(pvoxFname,' ',2);
nElec=size(pvoxCoordStr,1);
pvoxCoord=zeros(nElec,3);
for a=1:nElec,
    for b=1:3,
        pvoxCoord(a,b)=str2num(pvoxCoordStr{a,b});
    end
end
% Need to swap first two dimensions and flip 3rd to make coordinates
% compatible with vox2Seg.m
pvoxCoord=round(pvoxCoord+1);
pvoxCoord(:,[1 2])=pvoxCoord(:,[2 1]);
pvoxCoord(:,3)=sVol(3)-pvoxCoord(:,3);

% Import electrode labels
labelFname=fullfile(fsDir,subj,'elec_recon',sprintf('%s.electrodeNames',subj));
elecLabels=csv2Cell(labelFname,' ',2);

elecParc=cell(nElec,2);
hem=[];
for hemLoop=1:2,
    if hemLoop==1
        hem='L';
    else
        hem='R';
    end
    
    %% Are there any electrodes in this hemisphere?
    elecIdsThisHem=findStrInCell(hem,elecLabels(:,3));
    nElecThisHem=length(elecIdsThisHem);
    if nElecThisHem,
        %% READ SURFACE
        surfFname=fullfile(surfaceFolder,[lower(hem) 'h.pial']);
        [cort.vert, cort.tri]=read_surf(surfFname);
        nVertex=length(cort.vert);
        
        %% Get cortical parcellation
        switch upper(atlas)
            case 'DK'
                parcFname=fullfile(labelFolder,[lower(hem) 'h.aparc.annot']);
                [~, label, colortable]=read_annotation(parcFname);
                %[averts,label,colortable]=read_annotation(parcFname);
            case 'D'
                parcFname=fullfile(labelFolder,[lower(hem) 'h.aparc.a2009s.annot']);
                [~, label, colortable]=read_annotation(parcFname);
            case 'Y7'
                parcFname=fullfile(labelFolder,[lower(hem) 'h_Yeo2011_7Networks_N1000.mat']);
                load(parcFname);
            case 'Y17'
                parcFname=fullfile(labelFolder,[lower(hem) 'h_Yeo2011_17Networks_N1000.mat']);
                load(parcFname);
            otherwise
                error('Unrecognized value of atlas argument.')
        end
        
        for elecLoop=1:nElecThisHem,
            elecParc{elecIdsThisHem(elecLoop),1}=elecLabels{elecIdsThisHem(elecLoop),1};
            
            % Go through and set depth electrode assignments to depth:
            if elecLabels{elecIdsThisHem(elecLoop),2}=='D'
                elecParc{elecIdsThisHem(elecLoop),2}=vox2Seg(pvoxCoord(elecIdsThisHem(elecLoop),:),subj);
                %elecParc{elecIdsThisHem(elecLoop),2}='Depth'; % Could make
                %this optional if one only wants to use surface atlases
            else
                % Find closest vertex
                [~, minId]=min(sum( (repmat(pialCoord(elecIdsThisHem(elecLoop),:),nVertex,1)-cort.vert).^2,2 ));
                
                % Grab parcellation label for that vertex
                elecParc{elecIdsThisHem(elecLoop),2}=colortable.struct_names{find(colortable.table(:,5)==label(minId))};
            end
        end
        
    end
end

