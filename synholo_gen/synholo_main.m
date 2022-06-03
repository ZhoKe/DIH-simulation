% Synthetic hologram generation from 1) seeded particles and 2) a central 
% point scattered for PSF in deconvolution. 
% Written by Ke Zhou, 05/25/2022.

clear;clc;
% define particle field properties
Nf = 1; % numer of frames
Obj.enlargement = 2; % times that the holo being enlarged for accurate FFT
Obj.Supersmpl = 3; % supersampling rate of the obj and holo plane
Obj.size_obj = 256; % 256 descretization of the Obj, mesh refinement by a factor of "Supersmpl"
Obj.resolution = 5; % 5um/pixl; the SAME for both the obj and holo planes
Obj.dz = 5; % 5um/pixl; the SAME for both the obj and holo planes
Obj.z0 = 10e3; % 1e3 um, the distance from the hologram to the particle field center
Obj.Lz = 1280; % 1.6e3 um, depth of the measurement domain
Obj.n = 1; % refractive index
Obj.wavelength0 = 0.632; % um,  in the vacuum
Obj.rp = 4; %um, particle radius, note: diameter(2*rp) > resolution should be satisfied
Obj.Np = 1; % number of particles

xlim = [10, Obj.size_obj * Obj.resolution-10];
ylim = [10, Obj.size_obj * Obj.resolution-10];
zlim = [Obj.z0 - Obj.Lz/2 + 10, Obj.z0 + Obj.Lz/2-10];

%% particle seeding
position = cell(1,Nf);
for i = 1:Nf

[Xp, Yp, Zp] = P_seeding(xlim, ylim, zlim, Obj.Np);
[Zp, sortIdx] = sort(Zp);
position{i} = [(1:Obj.Np)', Xp(sortIdx), Yp(sortIdx), Zp];
% position{2} = [(1:Obj.Np)', Xp(sortIdx), Yp(sortIdx), Zp + 50];
end
%% synthetic hologram of particles
for i = 1:Nf
    disp(['holo simulation frame' num2str(i)])
    Obj.case = i;
    Obj.P_location = position{i};
    % hologram simulation
    [holo] = holo_gen(Obj);
    
end

%% synthetic hologram of a central point scatterer (PSF calculation)
Xp = (Obj.size_obj * Obj.resolution)/2 - Obj.resolution /2; % voxel coordinates of the point scatterer
Yp = Xp;
Zp = Obj.z0;

Obj.Np = 1; % number of particles
Obj.case = 0; % index = 0 for hologram of the single scatterer
Obj.rp = 4;  %um, particle radius, smaller than half voxel size
Obj.P_location = [1, Xp, Yp, Zp];

[holo_singlescatter] = holo_gen(Obj);

