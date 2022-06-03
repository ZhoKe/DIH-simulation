function [Np] = merge3D_all(Niteration,frame_i,volsiz,Z_depth, bwDrec_img, ALLbwDrec_img, obj_file,...
    L_merge_all, Dof_particle, par_list,show_particle_field_info)
% This function combine results from all iterations of IIPE with two methods.
% 1) combine recon particle fields, followed by particle extraction.
% 2) simply combine extracted particle centroids from all iterations.

%% 1) combine deconvolved reconstructions from all IIPE iterations
Nzscans = length(Z_depth);
for nz =1:Nzscans
    allbwDrec_img = double(imread(sprintf(bwDrec_img,frame_i,1,nz)));

    for iter =2:Niteration
        temp = double(imread(sprintf(bwDrec_img,frame_i,iter,nz)));
        allbwDrec_img = max (allbwDrec_img,temp);
    end
    imwrite(uint8(allbwDrec_img),sprintf(ALLbwDrec_img,frame_i,nz), 'tif', 'Compression', 'none' );

end

% 1.1) dilate and merge the binarized reconstructed optical field
BW = zeros(volsiz);

for nz = 1:volsiz(3)
    BW(:,:,nz) = imread(sprintf(ALLbwDrec_img,frame_i,nz));
end

BW = boolean(BW);
SE = strel('arbitrary',ones(L_merge_all,L_merge_all,L_merge_all));
BW3 = imdilate(BW,SE);
BW3 = bwlabeln(BW3);
BW3 = BW3.*BW;

% 1.2) translate pixles into voxel Corrdinate (x,y,z)
obj2 = regionprops(BW3,'PixelList');
[obj2.xyz] = obj2.PixelList;
obj2 = rmfield(obj2,'PixelList');

% 1.3) size filter to eliminate fake particles
index = zeros(1,length(obj2)); % preallocate a vector to store indexes of fake particles

for np = 1:length(obj2)
    Lp = max(obj2(np).xyz(:,3)) - min(obj2(np).xyz(:,3));

    if Lp < Dof_particle || isempty(obj2(np).xyz)
        index(np) = np;
    end

end
index = nonzeros(index); % extract nonzero indexes
obj2(index) = []; % delete fake particle indexes

% 1.4) merge close Centroids
Np = length(obj2);
xct = zeros(Np,1); xt1 = xct; xt2 = xct;
yct = zeros(Np,1); yt1 = yct; yt2 = yct;
zct = zeros(Np,1); zt1 = zct; zt2 = zct;

for np = 1:Np
    xct(np) = mean(obj2(np).xyz(:,1),1);
    yct(np) = mean(obj2(np).xyz(:,2),1);
    zct(np) = mean(obj2(np).xyz(:,3),1);

    xt1(np) = min(obj2(np).xyz(:,1));
    xt2(np) = max(obj2(np).xyz(:,1));
    yt1(np) = min(obj2(np).xyz(:,2));
    yt2(np) = max(obj2(np).xyz(:,2));
    zt1(np) = min(obj2(np).xyz(:,3));
    zt2(np) = max(obj2(np).xyz(:,3));

end

fprintf('After final length filter: %d objects\n', Np);

% 1.5) output particle locations
outfile = fopen(sprintf(par_list,frame_i),'w+');
fprintf(outfile,'TITLE = "3-D Particle Field"\r\n');
fprintf(outfile,'VARIABLES = " X ", " Y ", " Z ", " Xmin ", " Xmax ", " Ymin ", " Ymax ", " Zmin ", " Zmax "\r\n');
for np = 1:Np
    fprintf(outfile,'%f %f %f %f %f %f %f %f %f\r\n',xct(np),yct(np),zct(np), xt1(np),xt2(np), yt1(np),yt2(np), zt1(np),zt2(np) );
end
fclose(outfile);

%% 2) combine extracted particle centroids from all IIPE iterations
obj2_cmb = struct([]);
for i = 1:Niteration
    load(sprintf(obj_file, frame_i, i), 'obj2');
    obj2_cmb = cat(1, obj2, obj2_cmb);
end
obj2 = obj2_cmb;

save(sprintf(obj_file, frame_i, 0),'obj2','-v7.3');

%% 3D rendering of all particles obtained with Method 2
if show_particle_field_info ==1
    for np =1:length(obj2)
        if ~isempty(obj2(np).xyz)
            plot3(obj2(np).xyz(:,3),obj2(np).xyz(:,1),obj2(np).xyz(:,2),'.','color',[0.5 0.5 0.5]);
            hold on;
        end
    end
    plot3(zct,xct,yct,'blueo','MArkerFaceColor','blue')
    grid on;
    hold on;
    xlabel('y');
    ylabel('x');
    zlabel('z');
    title('blue: Geometrical centers');
end



end