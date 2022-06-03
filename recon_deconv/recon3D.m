% 3-D Reconstruction using a hologram
% Written by Ke Zhou, 05/25/2022

function [rec3] = recon3D(hologram, lambda, resolution, Z_depth)
ImgSiz = size(hologram);
rec3 = zeros(ImgSiz(1),ImgSiz(2),length(Z_depth));

%% Rayleight-Sommerfield kernel (in Fourier domain)
xn = 0:(ImgSiz(2)-1);
yn = 0:(ImgSiz(1)-1);
[X,Y] = meshgrid(xn,yn);
fx = (mod(X + ImgSiz(2)/2, ImgSiz(2)) - floor(ImgSiz(2)/2)) / ImgSiz(2);
fy = (mod(Y + ImgSiz(1)/2, ImgSiz(1)) - floor(ImgSiz(1)/2)) / ImgSiz(1);
f2 = fx.*fx + fy.*fy;
sqrt_input = 1 - f2*(lambda/resolution)^2;
sqrt_input(sqrt_input < 0) = 0;
H = 1i*(2*pi/lambda)*sqrt(sqrt_input);

%% Reconstruction slice-wise
FFT_holo = fft2(hologram);

for layer=1:length(Z_depth)

    % kernel formation
    Hz = exp(-Z_depth(layer).*H);

    % 3D recon with back-propagation
    rec = ifft2(FFT_holo.*Hz);
    rec3(:,:,layer) = rec(:,:);

end

end











