function [Nz] = sliceview(input, mytitle, mycolormap, mycaxis)

cmin = min(input(:)); % colormap min/max
cmax = max(input(:));
Nz = size(input,3); % number of layers along Z axis

figure();
for k = 1:Nz
subplot(2,ceil(Nz/2),k);
imagesc(squeeze(input(:,:,k))); axis image; caxis(mycaxis); colorbar; 
title([mytitle ' - layer ' num2str(k)]); xlabel('X');ylabel('Y'); colormap(mycolormap);
end

set(gcf,'position',[300 300 1200 500])

end