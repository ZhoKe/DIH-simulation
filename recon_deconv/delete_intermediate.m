function [flag] = delete_intermediate (output_path, flag)
% This function delete intermediate data generated in reconstruction.
% Use it after deconv-IIPE terminates for one hologram. Some intermediate results
% are used in deconv-IIPE and can not be deleted in computation.

if flag
fprintf('Deleting intermediate steps. This may take some time.\n');

delete([output_path, 'Dv_rec/*.tif']);
% rmdir([output_path, 'Dv_rec']);
% delete([output_path, 'Cmb_imgs/*']);
% rmdir([output_path, 'Cmb_imgs']);
delete([output_path, '3D_imgs/*.tif']);
% delete([output_path, '3D_imgs/*.mat']);

end

end