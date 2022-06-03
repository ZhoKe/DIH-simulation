function [clean_holo] = IIPE (iter, frame_i, obj_file, holo_img, reso, lambda, Z_depth, dz,...
    ImgSiz, C0, Niter, R_dial, Corast_recovery)

%% load data: particle position (obj_file) and hologram (holo_img)
p = load(sprintf(obj_file, frame_i,iter));
holo = (imread(sprintf(holo_img,frame_i,iter)));

%% local variables
Np_detected = length(p.obj2);

%% iterative particle extraction
for i = 1:Np_detected

    % voxel coordinates in a particle (vector)
    y = p.obj2(i).xyz(:,1);
    x = p.obj2(i).xyz(:,2);
    z = p.obj2(i).xyz(:,3);

    % particle removal
    [clean_holo] = particle_removel(holo, x, y, z, reso, lambda, Z_depth, dz, ImgSiz, C0,...
        Niter, R_dial, Corast_recovery);
    
    % store residual hologram for the next round of extraction
    holo = clean_holo;

end

%% save particle-extracted hologram
imwrite( uint8(abs(clean_holo)), sprintf( holo_img, frame_i,iter+1 ), 'tif', 'Compression', 'none' );

end