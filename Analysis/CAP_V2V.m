function Out=CAP_V2V(In,In_dim,In_mat,Out_dim,Out_mat)
% map voxels in the space of input volume to voxels in the space
% of output volumes
% Out_fn: uses only header info
%
% v1.0 Jonas Richiardi
% - initial release, based on code by Dimitri Van De Ville and Jonas
% Richiardi

% Out is filled with zeros and has the size of the output file
Out=zeros(Out_dim);

% generate all coordinates in output space
[x1,x2,x3]=ndgrid(1:Out_dim(1),1:Out_dim(2),1:Out_dim(3));
idx=1:numel(Out); % map all voxels

% take every voxel in the volume spanned by the output images,
% compute its real-world position in mm, then map input image

oobList=zeros(0,4); % list of out-of-bound input voxels

for iter=1:length(idx),
    oob=false;
    % recover world-space position of this voxel in mm from affine
    % transform matrix
    mm=Out_mat*[x1(idx(iter)) x2(idx(iter)) x3(idx(iter)) 1]';
    % convert this position into index of the closest structural voxel
    vx=round(In_mat\[mm(1) mm(2) mm(3) 1]'); 
    vx(vx<=0)=1;
    vxOri=vx;
    % remap out-of-bounds voxels to last  
    if vx(1)>In_dim(1), vx(1)=In_dim(1); oob=true; end 
    if vx(2)>In_dim(2), vx(2)=In_dim(2); oob=true; end
    if vx(3)>In_dim(3), vx(3)=In_dim(3); oob=true; end
    if (oob==true), oobList(end+1,:)=vxOri; end
    % idx(iter): current voxel
    Out(idx(iter))=In(vx(1),vx(2),vx(3));
    if any(Out(idx(iter))<0)  %Out
        warning('mapV2V:negativeVal',['Negative voxel values at ' num2str(iter)]);
    end
end;