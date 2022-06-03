% Optical Field Enhancement Using 3-D Deconvolution
% Written by: Santosh Kumar Sankar based on Dixon et.al 2011
% Flow Field Imaging Lab
% August 2015
function [Nz, avgI_deconv,stdI_deconv] = deconv3(holo_img, holo_psf_img, itern, Z_depth, frame_i,...
    lambda, reso, beta, nrm1, nrm2,avgI_deconv,stdI_deconv, deconv_on,...
    Drec_img, D_intensityfield, O_intensityfield, PSF_intensityfield,...
    OcmbXY_img, OcmbXZ_img, OcmbYZ_img, DcmbXY_img, DcmbXZ_img, DcmbYZ_img, flag_save)

%% Hologram Normalization
holo = double(imread(sprintf(holo_img,frame_i,itern)));
avg = mean(holo(:));
holo = holo./avg;
holo_n = holo - 1;

holo_psf = double(imread(holo_psf_img));
avg = mean(holo_psf(:));
holo_psf = holo_psf./avg;
holo_psf_n = holo_psf - 1;

%% 3-D Reconstructions (Optical Field & PSF):
[rec3] = recon3D(holo_n, lambda, reso, Z_depth);
[psf3] = recon3D(holo_psf_n ,lambda, reso, Z_depth);

%% Noise & DOF reduction using 3-D Deconvolution:
[Nx, Ny, Nz] = size(rec3);
rec3_int = rec3 .* conj(rec3) .* 10^( 2); % Scaled to get good SNR of the deconvoled field.
psf3_int = psf3 .* conj(psf3) .* 10^(-2); % Tuned empirically, but universal.

rec3_int_pad = padarray(rec3_int, [Nx/4, Ny/4, Nz/4], 0);
psf3_int_pad = padarray(psf3_int, [Nx/4, Ny/4, Nz/4], 0);

rec3_int_FT = fftn(rec3_int_pad);
psf3_int_FT = fftn(psf3_int_pad);

CTFi = conj(psf3_int_FT) ./ (psf3_int_FT .* conj(psf3_int_FT) + beta);
object_FT = rec3_int_FT .* CTFi;
object = abs(fftshift(ifftn(object_FT)));

object = object(0.75*Nx-Nx/2+1:0.75*Nx+Nx/2, 0.75*Ny-Ny/2+1:0.75*Ny+Ny/2, 0.75*Nz-Nz/2+1:0.75*Nz+Nz/2);

rec3_int_deconv = imtranslate(object, [-1, -1, -1]); % restore the lateral position after deconv

if ~deconv_on % cancel deconvolution
    rec3_int_deconv = rec3_int;
end
%% 1) save deconvolved plane images (have-to)
Nz = length(Z_depth);

if nrm1==1 % intensity normalization within the frame
    rec3_int_deconv  = (rec3_int_deconv - min(rec3_int_deconv(:)))./(max(rec3_int_deconv(:))-min(rec3_int_deconv(:)));
end

if nrm2==1 % intensity normalization across frames
    if itern ==1
        avgI_deconv = (mean(rec3_int_deconv(:)));
        stdI_deconv = (std(rec3_int_deconv(:)));
    elseif itern >1
        I_H2 = (mean(rec3_int_deconv(:)));
        S_H2 = (std(rec3_int_deconv(:)));
        rec3_int_deconv = avgI_deconv + (rec3_int_deconv - I_H2).*(stdI_deconv/S_H2);
    end
end

for i =1:Nz % save normalized deconvolved intensity field slice-wise as 8bit images
    CR = rec3_int_deconv(:,:,i);
    iCR = imcomplement(CR);
    iCR = iCR.*255;
    imwrite(uint8(iCR),sprintf(Drec_img, frame_i,itern,i), 'tif', 'Compression', 'none' );
end

%% 2) save min intensity maps in the X-Y plane (have-to)
rec3_int_deconv = imcomplement(rec3_int_deconv);
b = min(rec3_int_deconv,[],3);
imwrite(uint8(b*255),sprintf(DcmbXY_img,frame_i,itern));

%% 3) save results for diagnosis and plotting (optional)

if flag_save==1
    % Save optical fields
%     save(sprintf(D_intensityfield,frame_i,itern ),'rec3_int_deconv','-v7.3');
%     save(sprintf(O_intensityfield,frame_i,itern ),'rec3_int','-v7.3');
%     save(sprintf(PSF_intensityfield,frame_i,itern ),'psf3_int','-v7.3');

    % save min intensity maps for deconvolved fields in the XZ and YZ planes
    b = reshape(min(rec3_int_deconv,[],1),[size(rec3_int_deconv,2),size(rec3_int_deconv,3)]);
    imwrite(uint8(b*255),sprintf(DcmbXZ_img,frame_i,itern));

    b = reshape(min(rec3_int_deconv,[],2),[size(rec3_int_deconv,1),size(rec3_int_deconv,3)]);
    imwrite(uint8(b*255),sprintf(DcmbYZ_img,frame_i,itern));
    
    % save min intensity maps for original fields 
    rec3_int = (rec3_int - min(rec3_int(:)))./(max(rec3_int(:))-min(rec3_int(:)));
    rec3_int = imcomplement(rec3_int);
    b = min(rec3_int,[],3);
    imwrite(uint8(b*255),sprintf(OcmbXY_img,frame_i,itern));

    b = reshape(min(rec3_int,[],1),[size(rec3_int,2),size(rec3_int,3)]);
    imwrite(uint8(b*255),sprintf(OcmbXZ_img,frame_i,itern));

    b = reshape(min(rec3_int,[],2),[size(rec3_int,1),size(rec3_int,3)]);
    imwrite(uint8(b*255),sprintf(OcmbYZ_img,frame_i,itern));
end

end




