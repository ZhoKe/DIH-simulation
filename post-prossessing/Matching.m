
function [index_rec_ghost, index_rec_ambiguity, index_rec_match, index_gt_match, index_gt_ambiguity] = Matching(Nf, GT_position, Rec_position, Num_recparticle, Lz, Lr)

%% rec/gt matching
% index of three types of particles
index_rec_ghost = cell(Nf,1);
index_rec_ambiguity = cell(Nf,1); % one rec particle satisfies matching condition with several gt particles nearby
index_rec_match = cell(Nf,1);
% index of matched particles in ground truth
index_gt_match = cell(Nf,1);
index_gt_ambiguity = cell(Nf,1);

for i = 1:Nf
    for j = 1:Num_recparticle(i)  % for reconstructed particles
        
        % matching criterion calculation
        Z_dist = abs(Rec_position(i,j,3) - GT_position(i,:,3));
        r_dist = sqrt((Rec_position(i,j,1) - GT_position(i,:,1)).^2 + (Rec_position(i,j,2) - GT_position(i,:,2)).^2);
        indx_z = find(Z_dist<Lz);
        indx_r = find(r_dist<Lr);
        indx_match = intersect(indx_z, indx_r); % find indx of GT particles satisfying both criterion
        
        % matching flag and index storing
        % 1) for rec particles
        flag = length(indx_match);
        switch flag
            case 0 % ghost particles
                index_rec_ghost{i} = [index_rec_ghost{i}, j];
            case 1 % matched particles
                index_rec_match{i} = [index_rec_match{i}, j];
            otherwise % ambiguous particles
                index_rec_ambiguity{i} = [index_rec_ambiguity{i}, j];
        end
        
        % 2) for gt particles
        if flag>1 % if ambiguity occurs
            dist = zeros(length(indx_match),1);
            for k = 1:length(indx_match) % compare distances of the possible matches to the said rec particle
                dist(k) = sqrt((Rec_position(i,j,1) - GT_position(i,indx_match(k),1)).^2 ...
                    + (Rec_position(i,j,2) - GT_position(i,indx_match(k),2)).^2 ...
                    + (Rec_position(i,j,3) - GT_position(i,indx_match(k),3)).^2);
            end
            index_gt_ambiguity{i} = [index_gt_ambiguity{i}, indx_match(dist == min(dist))]; % find index of the closest matched gt particle
            clear dist
        else % if no ambiguity (flag == 1) or no matching (flag == 0)
            index_gt_match{i} = [index_gt_match{i}, indx_match];
        end   
        
    end
end

end