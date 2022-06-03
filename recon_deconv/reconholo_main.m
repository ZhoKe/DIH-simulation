clear;clc;

%% Parameters
% optical setup and flow domain
ImgSiz = [512 512]; % hologram pixel resolution
reso = 5; % [um] % lateral pixel resolution
dz = 5; % [um] % depth voxel resolution
z0 = 2e3; % [um] % object distance, from hologram to 3D domain center
Lz = 2560; % [um] % depth of 3D domain
lambda_vacu = 0.632; % [um] wavelength in vacuum
refractiveindex = 1.0; 
Dp = 2; % [pixels] estimated particle diameter

frameindex = [1:50];
IIPE_iterN = 2; % Number of IIPE iterations to be performed.
Dof_particle = 15; % [voxels] depth of the field. Min particle length

% flags
org = 0; % recon enabled or not
deconv_on = true;
show_particle_field_info =0;
flag_delete = 1;
flag_save_deconv = 1; % save intermediate results in the deconv step
flag_save_autothresh = 0; % save intermediate results in the auto-thresh step

% Infrequently changed parameters
beta = 10; % deconvolution relaxation parameter
R_dial = 2; % radius of dilation disk used in particle removal. Default 2.
C0 = 1.0; % scaling factor to apply to mask replacement intensity value. Default 1.
Niter = 4; % Number of inner iterations for particle removal. Default 4.
Contrast_recovery = 1; % Flag: contrast recovery in the removal step. Default 1.

% Derived parameters
nz = Lz/dz; % number of Z-layers in recon
lambda = lambda_vacu/refractiveindex;
Z_depth = z0 - Lz/2 + (1:1:nz) * dz;
volsiz = [ImgSiz, length(Z_depth)];
L_merge = 2*Dp;  % approximation for merger distance each iteration
L_merge_all = L_merge; % approximation for merger distance of all iterations
% Nzscans = length(Z_depth);
% layers = 1:length(Z_depth);

%% Create folders and paths for data saving/exchange
% setup I/O directory paths
mydir  = pwd;
idcs   = strfind(mydir,'\');
mydir(idcs) = '/'; % change the direction of slash to avoid errors in sprintf function
outterdir = mydir(1:idcs(end)-1);
path = [outterdir '/data'];

input_path = [path, '/syn_holo'];
input_hologram = [input_path, '/Org_Hologram_%05d.tif'];
output_path = [path, '/Output/'];

% create local folders
if ~isfolder([output_path, 'Holograms'])
    mkdir(output_path, 'Holograms');
end
if ~isfolder([output_path, 'Org_rec']) && org
    mkdir(output_path, 'Org_rec');
end
if ~isfolder([output_path, 'Dv_rec'])
    mkdir(output_path, 'Dv_rec');
end
if ~isfolder([output_path, 'Cmb_imgs'])
    mkdir(output_path, 'Cmb_imgs');
end
if ~isfolder([output_path, '3D_imgs'])
    mkdir(output_path, '3D_imgs');
end

% create local paths
% holograms
holo_img = [output_path, 'Holograms/Holo_%04d_iter%03d.tif'];
holo_psf_img = [output_path, 'Holograms/Holo_psf.tif'];
% 2D recon slices of of Original and Deconvolved recon
Orec_img = [output_path, 'Org_rec/Orec_%03d_%01d_%03d.tif'];
Drec_img = [output_path, 'Dv_rec/Drec_%03d_%01d_%03d.tif'];
% 3D recon distribution of O and D
O_intensityfield = [output_path, '3D_imgs/O_OF_%03d_%01d.mat'];
D_intensityfield = [output_path, '3D_imgs/D_OF_%03d_%01d.mat'];
PSF_intensityfield = [output_path, '3D_imgs/PSF_%03d_%01d.mat'];
% 2D min intensity maps of O and D
OcmbXY_img = [output_path, 'Cmb_imgs/OcmbXY_%03d_%01d.tif'];
OcmbXZ_img = [output_path, 'Cmb_imgs/OcmbXZ_%03d_%01d.tif'];
OcmbYZ_img = [output_path, 'Cmb_imgs/OcmbYZ_%03d_%01d.tif'];
DcmbXY_img = [output_path, 'Cmb_imgs/DcmbXY_%03d_%01d.tif'];
DcmbXZ_img = [output_path, 'Cmb_imgs/DcmbXZ_%03d_%01d.tif'];
DcmbYZ_img = [output_path, 'Cmb_imgs/DcmbYZ_%03d_%01d.tif'];
% locally-normalized 3D recon of D
NDcmbXY_img = [output_path, 'Cmb_imgs/NDcmbXY_%03d_%01d.tif'];
NDcmbYZ_img = [output_path, 'Cmb_imgs/NDcmbYZ_%03d_%01d.tif'];
NDcmbXZ_img = [output_path, 'Cmb_imgs/NDcmbXZ_%03d_%01d.tif'];
NDrec_img = [output_path, 'Dv_rec/NDrec_%03d_%01d_%03d.tif'];
% binarized 3D recon of D
bwDrec_img = [output_path, 'Dv_rec/BWDrec_%03d_%01d_%03d.tif'];
ALLbwDrec_img = [output_path, '3D_imgs/BW_%04d_%04d.tif'];
% extracted particle positions 
obj_file = [output_path, '3D_imgs/OBJs_%03d_%01d.mat']; % .mat file. For different iterations of IIPE
par_list = [output_path, '3D_imgs/ParticleCentroid_%03d.txt']; % .txt file. For the combined recons of D.

%% transfer synthetic holograms to the working directory
for i = 1:length(frameindex)
    copyfile(sprintf(input_hologram, 0), holo_psf_img);
    copyfile(sprintf(input_hologram, frameindex(i)), sprintf(holo_img, frameindex(i), 1));
end

% if estimate_DOF
%     if ~isdir([pathn, 'synth'])
%         mkdir(pathn, 'synth');
%     end
%     HoloSynthLz = [pathn, 'synth/Holo_Synth_Lz_001.tif'];
%     rec_imgSynth = [pathn, 'synth/rec_Synth_Lz_%03d_%03d.tif'];
%     LimPer = 0.75;
%
%     [Lz1,Lz,Lz2] = DOFestimator(HoloSynthLz,rec_imgSynth,ImgSiz, Dp, Lambda, Reso, Z_depth, LimPer);
%     Lz = round(Lz1/5);   % since Dv reduces DOF 5 times...
%
% end

%% START deconv + IIPE !
for frame_i = frameindex(1):frameindex(end)
    tic
    fprintf(['Reconstruction of frame' num2str(frame_i) '\n']);
    
    % preallocate variables
    avgI_deconv = 0; % mean of 3D deconvolution optical field for normalization
    stdI_deconv = 0; % std of 3D deconvolution optical field for normalization

    for iter = 1:IIPE_iterN

        if iter ==1
            itern =iter;
        end

        % STEP#1: 3D Reconstruction + Deconvolution:
        fprintf('Deconvolution...\n');
        [Nz_scans,avgI_deconv,stdI_deconv] = deconv3(holo_img, holo_psf_img, itern, Z_depth, frame_i, ...
                                       lambda, reso, beta, 1, 1, avgI_deconv, stdI_deconv, deconv_on,...
                                       Drec_img, D_intensityfield, O_intensityfield, PSF_intensityfield, ...
                                       OcmbXY_img, OcmbXZ_img, OcmbYZ_img, DcmbXY_img, DcmbXZ_img, DcmbYZ_img, flag_save_deconv);

        % STEP#2: 3D SNR Enhanceme + Thresholding:
        fprintf('Enhancement...\n');
        [thr_seg, thr_bg, BW] = auto_thresh(iter, Dp, Z_depth, frame_i, ImgSiz, DcmbXY_img,...
                                            Drec_img, NDcmbXY_img, NDcmbYZ_img, NDcmbXZ_img, NDrec_img, ...
                                            bwDrec_img, flag_save_autothresh);

        % STEP#3: Particle formation by merging 2D segments into 3D particle bodies
        fprintf('3D merging...\n');
        [no_more_object] = merge3D(iter,frame_i, BW, obj_file, L_merge, Dof_particle);

        % STEP#5: Iterative Inverse Particle Removal (IIPE)
        fprintf('Particle removal...\n');
        if no_more_object==0
            [Clean_Holo] = IIPE (iter, frame_i, obj_file, holo_img, reso, lambda, Z_depth, dz,...
                                 ImgSiz, C0, Niter, R_dial, Contrast_recovery);
            itern = iter+1;
        else
            itern =iter;
        end

        time = toc;
        fprintf(['Iteration ' num2str(iter) ' of frame ' num2str(frame_i) ' Done!  \n']);
    end

    % STEP#6: Merge the objects of all IIPEs:
    fprintf('Combining objects from all iterations...\n');
    [Np] = merge3D_all(IIPE_iterN, frame_i, volsiz, Z_depth, bwDrec_img, ALLbwDrec_img, obj_file,...
                                    L_merge_all, Dof_particle, par_list, show_particle_field_info);
    % STEP#7: Delete intermediate files
    fprintf('Deleting intermediate results...\n');
    [flag_delete] = delete_intermediate (output_path, flag_delete);

    toc
end




