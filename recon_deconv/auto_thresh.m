% Automatic Thresholding
% Written and Developed: Mostafa Toloui
% Flow Field Imaging Lab(Hong's group), ME-UMN
% April 2015

function [thr_seg, thr_bg, BW] = auto_thresh(iter, Dp, Z_depth, frame_i, ImgSiz , DcmbXY_img,...
                           Drec_img, NDcmbXY_img, NDcmbYZ_img, NDcmbXZ_img, NDrec_img, ...
                           bwDrec_img, flag_save)

%% local variable
IW1 = Dp*4; % interrogation window size for min intensity map scanning
IW1_step = IW1/2; % step length 

IW2 = 32; % block size of the min-max local normalization (75% of the particle DOF), subject to change 
IW2_step = IW2/2; % step length 
Nz = length(Z_depth); % number of depth layers

min_Imap = imread(sprintf(DcmbXY_img,frame_i,iter));

%% STEP1 : 2D IW scanning to calcualte thr0

c = zeros(size(min_Imap));
d = zeros(size(min_Imap));
e = zeros(size(min_Imap));

IW_index=1;
for i = 1:IW1_step:ImgSiz(1)-IW1+1
    for j = 1:IW1_step:ImgSiz(2)-IW1+1
        
        temp =  min_Imap(i:i+IW1-1,j:j+IW1-1);
        
        % min, max, avg and std for 2D IWs
        MN2(IW_index) = min(temp(:));
        MX2(IW_index) = max(temp(:));
        AVG2(IW_index) = mean2(temp(:));
        STD2(IW_index) = std2(temp(:));

        c(i:i+IW1-1,j:j+IW1-1) = (min_Imap(i:i+IW1-1,j:j+IW1-1) - MN2(IW_index))./(MX2(IW_index)-MN2(IW_index));
        d(i:i+IW1-1,j:j+IW1-1) = (min_Imap(i:i+IW1-1,j:j+IW1-1) - AVG2(IW_index))./(STD2(IW_index));
        e(i:i+IW1-1,j:j+IW1-1) = (STD2(IW_index));
        IW_index=IW_index+1;

    end
end

thr_bg = mean(double(MN2))-mean(double(STD2))-std((double(STD2)));

%% STEP 2: 3D thresholding using thr0
rec3 = zeros(ImgSiz(1),ImgSiz(2),Nz);
for iz = 1:Nz

    deconvslice = double(imread(sprintf(Drec_img,frame_i,iter,iz)));
    deconvslice(deconvslice>thr_bg) = 255;
    rec3(:,:,iz) = deconvslice;

end

% rec3(:,ImgSiz(1),:) = max(rec3(:));
% rec3(ImgSiz(1),:,:) = max(rec3(:));

%% STEP 3: 3D local intensity normalization
IW_index=1;
for k = 1:IW2_step:Nz-IW2+1
    for i = 1:IW2_step:ImgSiz(1)-IW2+1
        for j = 1:IW2_step:ImgSiz(2)-IW2+1
            
            blk =  rec3(i:i+IW2-1,j:j+IW2-1,k:k+IW2-1);

            MN3(IW_index) = (min(blk(:)));
            MX3(IW_index) = (max(blk(:)));
            ST3(IW_index) = (std(blk(:)));
            AVG3(IW_index) = (mean(blk(:)));

            if (MN3(IW_index) <= thr_bg)
                rec3_norm(i:i+IW2-1,j:j+IW2-1,k:k+IW2-1) = (blk - MN3(IW_index))./(MX3(IW_index)-MN3(IW_index));  
            else
                rec3_norm(i:i+IW2-1,j:j+IW2-1,k:k+IW2-1) = 1;
            end

            IW_index = IW_index+1;

        end
    end
end
rec3_norm = rec3_norm * 255;

%% Automatic thresholding 
bb = reshape(min(rec3_norm,[],3), ImgSiz(1)*ImgSiz(2),1);

To = mean(bb);
[counts,centers] = hist(bb, max(bb));
IW_index=1;

for i = round(To):-1:min(bb)+1
    if counts(i)~=0
        thr_list(IW_index) = centers(i);
        IW_index=IW_index+1;
    end
end

thr_seg = max(thr_list); % threshold calculation

BW = zeros(size(rec3_norm)); % BW: binarized 3D domain
BW(rec3_norm < thr_seg) = 1; % thresholding the locally-normalized 3D recon for particle segmentation

%% save intermediate results
for iz = 1:Nz % (have-to)

    slice_BW = BW(:,:,iz);
    imwrite( (slice_BW) , (sprintf(bwDrec_img,frame_i,iter,iz)), 'tif', 'Compression', 'none' );

end

if (flag_save ==1) % optional
    clear temp;
    temp = min(rec3_norm,[],3);
    imwrite( uint8(temp) , (sprintf(NDcmbXY_img,frame_i,iter)), 'tif', 'Compression', 'none' );
    temp = reshape(min(rec3_norm,[],1),[size(rec3_norm,2),size(rec3_norm,3)]);
    imwrite( uint8(temp) , (sprintf(NDcmbXZ_img,frame_i,iter)), 'tif', 'Compression', 'none' );
    temp = reshape(min(rec3_norm,[],2),[size(rec3_norm,1),size(rec3_norm,3)]);
    imwrite( uint8(temp) , (sprintf(NDcmbYZ_img,frame_i,iter)), 'tif', 'Compression', 'none' );

    % BW and rec3_norm
    for iz = 1:Nz

        slice_rec3_norm = rec3_norm(:,:,iz);

        imwrite( uint8(slice_rec3_norm) , (sprintf(NDrec_img,frame_i,iter,iz)), 'tif', 'Compression', 'none' )

    end

end

end
