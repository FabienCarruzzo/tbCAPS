%% Converts CAPs (matlab matrix) into NIFTI files
% CAP must have size n_CAPs x n_voxels
function [] = CAPToNIFTI(CAP,mask,brain_info,savedir,savename)

    % Number of CAPs
    n_CAPs = size(CAP,1);

    % Voxel size
    voxel_size = diag(brain_info.mat);
    voxel_size = voxel_size(1:end-1)';
    
    voxel_shift = brain_info.mat(:,4);
    voxel_shift = voxel_shift(1:end-1)';
    
    % Converts each CAP into a 3D volume
    for i = 1:n_CAPs
        
        tmp = CAP(i,:);
        V = zeros(brain_info.dim);
        V(mask) = tmp;

        tmp_NIFTI = make_nii(V,voxel_size,-voxel_shift./voxel_size);
        
        tmp_NIFTI.hdr.dime.datatype=64;
        tmp_NIFTI.hdr.dime.bitpix=64;
        
        %y_Write(V,brain_info,fullfile(savedir,[savename,'_CAP',num2str(i),'.nii']));
        save_nii(tmp_NIFTI,fullfile(savedir,[savename,'_CAP',num2str(i),'.nii']));
    end
end