clc;clear;

%% I/O directory paths
mydir  = pwd;
idcs   = strfind(mydir,'\');
mydir(idcs) = '/'; % change the direction of slash to avoid errors in sprintf function
outterdir = mydir(1:idcs(end)-1);
path = [outterdir '/data'];

GT = [path '/syn_holo/GT_position_%03d.mat'];
Rec = [path '/Output/3D_imgs/OBJs_%03d_0.mat'];

%% parameters
frameindex = [1:50];
Nf = length(frameindex); % number of frames
resolution = 5; % [um]
z0 = 2e3; % [um] the distance from the hologram to the particle field center
Lz = 2560; % [um] depth of the measurement domain
z_start = (z0 - Lz/2)/resolution; % [voxel] 
Lz = 16; % searching radius along Z axis
Lr = 1; % searching radius in the lateral plane

% define array to store particle positions
GT_position = []; % dim1 = frame; dim2 = particle; dim3 = coordinates
Rec_position = []; % dim1 = frame; dim2 = particle; dim3 = coordinates

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

%% Matching
[index_rec_ghost, index_rec_ambiguity, index_rec_match, index_gt_match, index_gt_ambiguity] = Matching(Nf, GT_position, Rec_position, Num_recparticle, Lz, Lr);

%% Meas. PDF
for i = 1:Nf
dx = Rec_position(i,index_rec_match{i},1) - GT_position(i,index_gt_match{i},1);
dy = Rec_position(i,index_rec_match{i},2) - GT_position(i,index_gt_match{i},2);
dz = Rec_position(i,index_rec_match{i},3) - GT_position(i,index_gt_match{i},3);

meanx(i) = mean(dx); stdx(i) = std(dx);
meany(i) = mean(dy); stdy(i) = std(dy);
meanz(i) = mean(dz); stdz(i) = std(dz);

figure(i);
subplot(2,2,1)
histogram(dx, 40); xlim([-1 1]); ylim([0 200]); xlabel('dx (vx)'); ylabel('Counts'); title('Meas. PDF - X'); % ylim([0 50]);
subplot(2,2,2)
histogram(dy, 40); xlim([-1 1]); ylim([0 200]); xlabel('dy (vx)'); ylabel('Counts'); title('Meas. PDF - Y'); % ylim([0 50]);
subplot(2,2,3)
histogram(dz, 40); xlim([-20 20]); ylim([0 200]); xlabel('dz (vx)'); ylabel('Counts'); title('Meas. PDF - Z'); % ylim([0 80]);
end

%% Summary of 3 types of rec particles (ghost, matched, ambg)
disp(['Num of total rec particles   = ' num2str(Num_recparticle)]);

num_particle_matched= cellfun(@length,index_rec_match);
disp(['Num of rec matched particles = ' num2str(num_particle_matched')]);

num_particle_ghost= cellfun(@length,index_rec_ghost);
disp(['Num of rec ghost particles   = ' num2str(num_particle_ghost')]);

num_particle_ambg= cellfun(@length,index_rec_ambiguity);
disp(['Num of rec ambg particles    = ' num2str(num_particle_ambg')]);

disp(['Stdx   = ' num2str(stdx) '    [vx]']);
disp(['Stdy   = ' num2str(stdy) '    [vx]']);
disp(['Stdz   = ' num2str(stdz) '    [vx]']);
disp(['Meanx  = ' num2str(meanx) '    [vx]']);
disp(['Meany  = ' num2str(meany) '    [vx]']);
disp(['Meanz  = ' num2str(meanz) '    [vx]']);

%% plot rec/gt data
for i = 1:Nf
    figure(i);
    subplot(2,2,4)
    plot3(GT_position(i,[index_gt_match{i} index_gt_ambiguity{i}],3),...
          GT_position(i,[index_gt_match{i} index_gt_ambiguity{i}],1),...
          GT_position(i,[index_gt_match{i} index_gt_ambiguity{i}],2),'ro','MArkerFaceColor','r', 'Markersize', 2);
    hold on;
    
%     plot3(GT_position(i,index_gt_ambiguity{i},3),GT_position(i,index_gt_ambiguity{i},1),GT_position(i,index_gt_ambiguity{i},2),'redo','MArkerFaceColor','r');
%     hold on;
    
    plot3(Rec_position(i,index_rec_match{i},3),Rec_position(i,index_rec_match{i},1),Rec_position(i,index_rec_match{i},2),'bo','MArkerFaceColor','b', 'Markersize', 2);
    hold on;
    
    plot3(Rec_position(i,index_rec_ghost{i},3),Rec_position(i,index_rec_ghost{i},1),Rec_position(i,index_rec_ghost{i},2),'ko','MArkerFaceColor','k', 'Markersize', 2);
    hold on;
    
    plot3(Rec_position(i,index_rec_ambiguity{i},3),Rec_position(i,index_rec_ambiguity{i},1),Rec_position(i,index_rec_ambiguity{i},2),'go','MArkerFaceColor','green', 'Markersize', 2);
    
    grid on; axis image;
    xlabel('z (voxel)'); xlim([-10 522]);
    ylabel('x (voxel)'); ylim([-10 522]);
    zlabel('y (voxel)'); zlim([-10 522]);
    title(['Frame' num2str(i)]);
    legend('GT-matched', 'Rec-matched', 'Rec-ghost', 'Rec-ambg')

end