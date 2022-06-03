function [clean_holo] = particle_removel (holo, x, y, z, reso, lambda, Z_depth, dz, ImgSiz, C0,... 
                                          Niter, R_dial, Contrast_recovery)

z = mean(z); % in-focus Z location, or mass center of the recon particle. Unit: [voxel].
L = (z - 1)*dz+Z_depth(1); % in-focus object distance 

%% generate mask of the in-focus particle cross section (dilated by R_dial)
mask = zeros(ImgSiz);

for i=1:length(x)
    mask (x(i),y(i)) = 1; % Mask equals to unit for the cross sectional area with a particle
end

SE = strel('disk', R_dial); % dilation element
mask = imdilate(mask, SE);

%% in-focus optical field reconstruction (backward propagation)
reconParam=struct( 'imgHolo',(holo) , 'wavelength', lambda, 'res', reso, ...
    'dist', -L, 'filtType','RS' , 'outType', 'complex'  );
[Ap] = digiReconEx03( reconParam ); % cross section of a reconstruction at in-focus distance of one particle

% figure(); imagesc(holo); axis image; colormap gray; title('original holo');
%% iterative particle removal
for i=1:Niter

    % 1) In-focus Optical field mean estimation
    I_bg = mean(Ap(:));

    % 2) fill the bg intensity into the in-focus particle cross section
    CF_Ap = ((Ap).*(1-mask))+ ((C0*I_bg).*mask);

%     figure(); imagesc(abs(Ap)); axis image; colormap gray; title('Ap - before extraction');
%     figure(); imagesc(abs(CF_Ap)); axis image; colormap gray; title('Ap - masked');

    % 3) clean synthetic hologram generation (forward propagation)
    reconParam=struct( 'imgHolo',(CF_Ap) , 'wavelength', lambda, 'res', reso, ...
        'dist', L, 'filtType','RS' , 'outType', 'complex'  );
    [clean_holo] = digiReconEx03( reconParam );
    clean_holo=abs(clean_holo);

    % 4) contrast recovery using the previous hologram as a reference
    if ( Contrast_recovery == 1)
        [clean_holo] = Contrast_recov (holo, clean_holo);
    end
%       figure(); imagesc(clean_holo); axis image; colormap gray; title('clean holo');

    % 5) in-focus optical field reconstruction (backward propagation)
    reconParam=struct( 'imgHolo',(clean_holo) , 'wavelength', lambda, 'res', reso, ...
        'dist', -L, 'filtType','RS' , 'outType', 'complex'  );
    [Ap] = digiReconEx03( reconParam );

%     figure(); imagesc(abs(Ap)); axis image; colormap gray; title('Ap - after extraction');

end

end


