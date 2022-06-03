function [hologram] = holo_gen(obj)
% Synthetic hologram simulation using the Rayleigh-Summerfeld diffraction
% kernel (namely the Angular Spectrum Method).
% Written by Ke Zhou, 05/25/2022.

% imported parameters
num_particle = obj.Np;
caseindex = obj.case;
enlargement = obj.enlargement; % times that the holo being enlarged for accurate FFT
Supersmpl = obj.Supersmpl; % supersampling rate of the obj and holo plane
size_obj = obj.size_obj * Supersmpl; % descretization of the Obj, mesh refinement by a factor of "Supersmpl"
resolution = obj.resolution/Supersmpl; % [um/pixl]; the SAME for both the obj and holo planes
z0 = obj.z0; % [um], the distance from the hologram to the particle field center
Lz = obj.Lz; % [um], depth of the measurement domain
n = obj.n; % refractive index
wavelength0 = obj.wavelength0; % [um],  in the vacuum
radi = obj.rp/resolution; % [voxel], particle radi
P_location = obj.P_location;

% derived parameters
wavelength = wavelength0/n;
size_holo = enlargement*size_obj;

% local variables
FFT_holo = zeros(size_holo, size_holo);
FFT_ref = zeros(size_holo, size_holo);
hologram = zeros(obj.size_obj, obj.size_obj);
ref_hologram = zeros(obj.size_obj, obj.size_obj);

% output directive
mydir  = pwd;
idcs   = strfind(mydir,'\');
outterdir = mydir(1:idcs(end)-1);
out_pathn = [outterdir '\data\syn_holo'];

if ~isfolder(out_pathn)
    mkdir(out_pathn);
end
holo_img = '/Org_Hologram_%05d.tif';
holo_ref_img = '/Org_Hologram_%05d_ref.tif';
GT_plocation = '/GT_position_%03d.mat';

%% Rayleigh-Sommerfeld kernel in F domain
% Constant frequencies. write it this way to avoid fftshifts
xn = 0:(size_holo-1);
yn = 0:(size_holo-1);
[X,Y] = meshgrid(xn,yn);
fx = (mod(X + size_holo/2, size_holo) - floor(size_holo/2)) / size_holo;
fy = (mod(Y + size_holo/2, size_holo) - floor(size_holo/2)) / size_holo;
f2 = fx.*fx + fy.*fy;

sqrt_input = 1 - f2*(wavelength/resolution)^2;
sqrt_input(sqrt_input < 0) = 0;

H = 1i*(2*pi/wavelength)*sqrt(sqrt_input);

%% Convolution of objects with the kernel to simulate diffraction wave propagation
for n = 1:num_particle
    disp(['particle' num2str(n)]);

    xp = P_location(n, 2)/resolution;
    yp = P_location(n, 3)/resolution;
    zp = P_location(n, 4);

    % generate the particle mask function (0 for particle blockage, 1 for vaccum transmittance)
    [x_grid,y_grid] = meshgrid(1:size_obj, 1:size_obj);
    x_grid = x_grid -0.5; % offset coordinates to the center of a pixel
    y_grid = y_grid -0.5;

    P_field = (((y_grid - yp).^2+(x_grid - xp).^2) ./ (radi)^2)<= 1;
    P_field = padarray(P_field,[(size_holo-size_obj)/2  (size_holo-size_obj)/2], 'both'); % zero padding to match with the fftsize
    P_field = 1- P_field; %

    % generate reference mask function (unity matrix for vaccum)
    %R_field = ones(size(P_field,1), size(P_field,1)); % reference field

    % calculate phase
    z = zp; %-z
    phase = exp(1i*2*pi*(z0+Lz/2-z)/wavelength);

    % calculate FFT of the obj 
    FFT_obj = fft2(P_field*phase);
    %FFT_obj_ref = fft2(R_field*phase);

    % sum all propagated wave in the Fourier domain
    Hz = exp(z*H);
    FFT_holo = FFT_holo + (FFT_obj .* Hz); % negative/positive for dark/bright objects
    %FFT_ref = FFT_ref + (FFT_obj_ref .* Hz);
    
end

holo = ifft2(FFT_holo);
%ref = ifft2(FFT_ref);

%% crop the central reion if hologram was enlarged by > 1
if enlargement ~=1
    holo = holo(size_holo/2-size_obj/2+1:size_holo/2+size_obj/2,...
        size_holo/2-size_obj/2+1:size_holo/2+size_obj/2);
    
    %ref = ref(size_holo/2-size_obj/2+1:size_holo/2+size_obj/2,...
        %size_holo/2-size_obj/2+1:size_holo/2+size_obj/2);
end

%% Actual recorded hologram intensity (high-resolution)
holo_intensity = holo .* conj(holo);
%ref_intensity = ref .* conj(ref);

%% Downsized hologram (low-resolution)
for i = 1:obj.size_obj
    for j = 1:obj.size_obj
        hologram(i, j) = mean(holo_intensity((i-1)*Supersmpl+1:i*Supersmpl,...
                                             (j-1)*Supersmpl+1:j*Supersmpl),'all');
                                         
        %ref_hologram(i, j) = mean(ref_intensity((i-1)*Supersmpl+1:i*Supersmpl,...
                                             %(j-1)*Supersmpl+1:j*Supersmpl),'all');                                 
    end
end

%% normalize and output
hologram_norm = (hologram - min(hologram(:)))./(max(hologram(:)) - min(hologram(:)));
%ref_hologram_norm = (ref_hologram - min(hologram(:)))./(max(hologram(:)) - min(hologram(:)));

imwrite(uint8(hologram_norm*2^8), [out_pathn, sprintf(holo_img, caseindex)], 'tif', 'Compression', 'none');
%imwrite(uint8(ref_hologram_norm*2^8),[out_pathn, sprintf(holo_ref_img, caseindex)], 'tif', 'Compression', 'none');

save([out_pathn, sprintf(GT_plocation, caseindex)], 'P_location')

end

