clc;clear;

%% parameters

% DIH parameters
frameindex = [1:50];
Nf = length(frameindex); % number of frames
resolution = 5; % [um] resolution of discretization in both lateral and depth direction (isotropic)
z0 = 2e3; % [um] the distance from the hologram to the particle field center
Lz = 2560; % [um] depth of the measurement domain
Nx = 512; Ny = 512; Nz = 512; % discretization of the domain
z_start = (z0 - Lz/2)/resolution; % [voxel] starting point of the 3D domain
Lz = 16; % [voxel] searching radius along Z axis
Lr = 1; % [voxel] searching radius in the lateral plane

% sub-domain parameters
blocksize = round(Nx/6); % size of each sub-domain in voxels
viewon = 1; % view on/off

% define array to store particle positions
GT_position = []; % dim1 = frame; dim2 = particle; dim3 = coordinates
Rec_position = []; % dim1 = frame; dim2 = particle; dim3 = coordinates

%% I/O directory paths
mydir  = pwd;
idcs   = strfind(mydir,'\');
mydir(idcs) = '/'; % change the direction of slash to avoid errors in sprintf function
outterdir = mydir(1:idcs(end)-1);
path = [outterdir '/data'];

GT = [path '/syn_holo/GT_position_%03d.mat'];
Rec = [path '/Output/3D_imgs/OBJs_%03d_0.mat'];

%% Read rec/gt data
for i = 1:Nf
    load(sprintf( GT, frameindex(i)));
    GT_position(i,:,:) = P_location(:,2:4); % column 2 to 4 are x, y, z coordinates
    
    load(sprintf( Rec, frameindex(i)));
    for j = 1:length(obj2)
        Rec_position(i,j,:) = mean(obj2(j).xyz(:,:),1);
        Num_recparticle(i) =  length(obj2);
    end
end

Rec_position(:,:,1:2) = Rec_position(:,:,1:2) - 0.5; % offset coordinates to the center of voxels in the lateral plane

GT_position(:,:,1:2) = GT_position(:,:,1:2)/resolution;
GT_position(:,:,3) = GT_position(:,:,3)/resolution - z_start; 

%% Matching (per frame)
[index_rec_ghost, index_rec_ambiguity, index_rec_match, index_gt_match, index_gt_ambiguity] = Matching(Nf, GT_position, Rec_position, Num_recparticle, Lz, Lr);


%% particle screening and pdf calculation in the box region (for all frames)

% partition the 3D voxel domain

[X, Y, Z] = meshgrid(0:blocksize:Nx, 0:blocksize:Ny, 0:blocksize:Nz);

% preallocate std and mean matrices
STDX = zeros(size(X)-1); STDY = zeros(size(X)-1); STDZ = zeros(size(X)-1);
MEANX = zeros(size(X)-1); MEANY = zeros(size(X)-1); MEANZ = zeros(size(X)-1);

% Start block-wise operation
for k = 1:(size(Z,3)-1)
    for i = 1:(size(Y,1)-1)
        for j = 1:(size(X,2)-1)

            % vertice coordinates of the box
            X_box = [X(i,j,k), X(i,j+1,k)];
            Y_box = [Y(i,j,k), Y(i+1,j,k)];
            Z_box = [Z(i,j,k), Z(i,j,k+1)];

            % screen out particles inside the box
            [GT_boxed, Rec_boxed] = box(X_box, Y_box, Z_box, index_rec_match,...
                                        index_gt_match, GT_position, Rec_position, Nf);
            
            % view pdf and statistics
            [stdx, stdy, stdz, meanx, meany, meanz] = pdfview(GT_boxed, Rec_boxed, viewon); 

            % save statistics
            STDX(i,j,k) = stdx;
            STDY(i,j,k) = stdy;
            STDZ(i,j,k) = stdz;
            MEANX(i,j,k) = meanx;
            MEANY(i,j,k) = meany;
            MEANZ(i,j,k) = meanz;

        end
    end
end

%% slice plots of layers along the depth direction
[Nz] = sliceview(STDX, 'StdX', 'parula', [0.18 0.28]);
[Nz] = sliceview(STDY, 'StdY', 'parula', [0.18 0.28]);
[Nz] = sliceview(STDZ, 'StdZ', 'parula', [3 5.5]);
[Nz] = sliceview(MEANX, 'MeanX', 'jet', [-0.04 0.04]);
[Nz] = sliceview(MEANY, 'MeanY', 'jet', [-0.04 0.04]);
[Nz] = sliceview(MEANZ, 'MeanZ', 'jet', [-3.2 3.2]);