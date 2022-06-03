function [stdx, stdy, stdz, meanx, meany, meanz] = pdfview(GT, Rec, viewon)
% This function calculate PDF of particle localization errors.
% Input: 
% GT and Rec, groundtruth and reconstructed particle positions in voxel.
% viewon: logical variable to control the plot.
% Output: std and mean of the position errors (in voxel) along x, y and z directions.
%% Meas. PDF
dx = Rec(:,1) - GT(:,1);
dy = Rec(:,2) - GT(:,2);
dz = Rec(:,3) - GT(:,3);

meanx = mean(dx); stdx = std(dx);
meany = mean(dy); stdy = std(dy);
meanz = mean(dz); stdz = std(dz);




%% Summary of stats

Np_matched= length(GT);
disp(['Num of rec matched particles = ' num2str(Np_matched')]);

disp(['Stdx   = ' num2str(stdx) '    [vx]']);
disp(['Stdy   = ' num2str(stdy) '    [vx]']);
disp(['Stdz   = ' num2str(stdz) '    [vx]']);
disp(['Meanx  = ' num2str(meanx) '    [vx]']);
disp(['Meany  = ' num2str(meany) '    [vx]']);
disp(['Meanz  = ' num2str(meanz) '    [vx]']);



%% plot pdf and 3D distributions of rec/gt particles
if viewon
    figure(1);
    % three pdfs
    subplot(2,2,1)
    histogram(dx, 20); xlim([-1 1]); %ylim([0 200]);
    xlabel('dx (vx)'); ylabel('Counts'); title('Meas. PDF - X'); % ylim([0 50]);
    subplot(2,2,2)
    histogram(dy, 20); xlim([-1 1]); %ylim([0 200]);
    xlabel('dy (vx)'); ylabel('Counts'); title('Meas. PDF - Y'); % ylim([0 50]);
    subplot(2,2,3)
    histogram(dz, 20); xlim([-20 20]); %ylim([0 200]);
    xlabel('dz (vx)'); ylabel('Counts'); title('Meas. PDF - Z'); % ylim([0 80]);

    % particle distributions
    subplot(2,2,4)
    plot3(GT(:,3), GT(:,1), GT(:,2),'ro','MArkerFaceColor','r', 'Markersize', 1);
    hold on;
    plot3(Rec(:,3), Rec(:,1), Rec(:,2),'bo','MArkerFaceColor','b', 'Markersize', 1);

    grid on; axis image;
    xlabel('z (voxel)'); xlim([-10 522]);
    ylabel('x (voxel)'); ylim([-10 522]);
    zlabel('y (voxel)'); zlim([-10 522]);
    title(['Particle distribution in the box region']);
    legend('GT-matched', 'Rec-matched')
end

end