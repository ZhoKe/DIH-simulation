% 3D Particle Formation
% Written and Developed: Mostafa Toloui
% Flow Field Imaging Lab(Hong's group), ME-UMN
% April 2015

function [no_more_object] = merge3D(iter,frame_i, BW, obj_file, h, Dof_particle)

h1 = waitbar(0);
%% dilate and merge the binarized reconstructed optical field
BW = boolean(BW);
SE = strel('arbitrary',ones(h,h,h));
BW3 = imdilate(BW,SE);
BW3 = bwlabeln(BW3);
BW3 = BW3.*BW;

%% translate pixles into voxel Corrdinate (x,y,z)
obj2 = regionprops(BW3,'PixelList');
[obj2.xyz] = obj2.PixelList;
obj2 = rmfield(obj2,'PixelList');

%% size filter to eliminate fake particles
index = zeros(1,length(obj2)); % preallocate a vector to store indexes of fake particles

for i = 1:length(obj2)
    Lp = max(obj2(i).xyz(:,3)) - min(obj2(i).xyz(:,3));

    if Lp < Dof_particle || isempty(obj2(i).xyz)
        index(i) = i;
    end

end
index = nonzeros(index); % extract nonzero indexes
obj2(index) = []; % delete fake particle indexes

%% Save
if (~isempty(obj2))
    save(sprintf(obj_file, frame_i, iter),'obj2','-v7.3');
    no_more_object = 0;
else
    no_more_object = 1;
end

close(h1);


end