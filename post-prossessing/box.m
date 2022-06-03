function [GT_boxed, Rec_boxed] = box(X, Y, Z, index_rec_match, index_gt_match, GT_position, Rec_position, Nf)
% This function screen out the gt and rec particles inside a predefined box
% region, where the particle localization pdf can be calculated.
% Input:
% X, Y, Z: min/max coordinates of the box region along three directions in voxels. 
% Dimension = [1,2]. Other inputs have been explained elsewhere.
% Output:
% GT_boxed, Rec_boxed: groundtruth and reconstructed particle positions
% within the box region in voxels, for all the frames. Dimension = [X, 3], X is the
% number of particles detected in this box region for all frames.

GT_total = [];
Rec_total = [];

% combine coordinates of matched gt/rec particles
for i = 1:Nf
GT_total = [GT_total; squeeze(GT_position(i,index_gt_match{i},:))];
Rec_total = [Rec_total; squeeze(Rec_position(i,index_rec_match{i},:))];

end

% screen out idx of gt particles inside the box
GT_total_x = GT_total(:,1);
GT_total_y = GT_total(:,2);
GT_total_z = GT_total(:,3);

idx_x = find(GT_total_x > X(1) & GT_total_x < X(2));
idx_y = find(GT_total_y > Y(1) & GT_total_y < Y(2));
idx_z = find(GT_total_z > Z(1) & GT_total_z < Z(2));

idx = intersect(idx_x, idx_y);
idx = intersect(idx_z, idx);

% generate gt/rec particle coordinates
GT_boxed = GT_total(idx,:);
Rec_boxed = Rec_total(idx,:);

end