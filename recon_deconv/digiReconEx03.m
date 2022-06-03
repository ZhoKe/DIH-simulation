function [reconImages] = digiReconEx03( reconParam )
% function for generating digital hologram using matlab
%   reconParam - struct( 'imgHolo', [], 'wavelength', [], 'res', [], 'dist' [], 'filtType',[] );
%       again. [col, row], currently we only dealing with sensor with aspec ratio of 1:1
%       filtType: [COS(cosine) | RS (Rayleigh-Sommerfield) | KF (Kirchoff-Fersnel)];
%       outType: [ Intensity | complex ]
% (1) creating memory storage
clear j;
holoImgSiz = size( reconParam.imgHolo );

colImg = 2^ceil( log(holoImgSiz(2))/log(2) );       %% enlarged data set
rowImg = 2^ceil( log(holoImgSiz(1))/log(2) );

[digiHoloTemplKL] = createHoloTempl( [holoImgSiz(2), holoImgSiz(1)] );

switch( reconParam.outType )
    case 'Intensity'
        reconImages = uint8( zeros([holoImgSiz(1), holoImgSiz(2), length( reconParam.dist )]) );
    case 'complex'
        reconImages = complex( zeros([holoImgSiz(1), holoImgSiz(2), length( reconParam.dist )]), ...
            zeros([holoImgSiz(1), holoImgSiz(2), length( reconParam.dist )]) );
    otherwise
end
phaseConst = reconParam.wavelength*pi/(reconParam.res(1))^2;

dhImage = zeros( [rowImg, colImg] );        % create computing template
dhImage( 1:holoImgSiz(1), 1:holoImgSiz(2) ) = im2double( reconParam.imgHolo );  % load current data

dhImgFFT = fft2( dhImage );                 % create fft2 for current image

for i=1:length( reconParam.dist )
    switch( reconParam.filtType )
        case 'COS'     %% Cosine kernel
            dhfilter = cos(reconParam.dist(i)*phaseConst*digiHoloTemplKL);
            dhfilter = complex( dhfilter, zeros(size(dhfilter)) );

        case 'RS'       %% Rayleigh-Sommerfield formula
%             disp( 'RS' );
            dhfilter0 = 1-(reconParam.wavelength/reconParam.res(1))^2*digiHoloTemplKL;
			idx =  dhfilter0 <= 0 ;
			dhfilter0(idx) = 0;
            dhfilter1 = sqrt( dhfilter0 );
            dhfilter = exp(-1i*2*pi*reconParam.dist(i)/reconParam.wavelength*dhfilter1);
            clear dhfilter1;

        case 'RSG'       %% Rayleigh-Sommerfield formula - Grier Version
            disp( 'RSG' );
            dhfilter0 = 1-((reconParam.wavelength/reconParam.res(1))^2)*digiHoloTemplKL;
			idx =  dhfilter0 <= 0 ;
			dhfilter0(idx) = 0;
            dhfilter1 = sqrt( dhfilter0 );
            dhfilter = exp(1i*2*pi*reconParam.dist(i)/reconParam.wavelength*(dhfilter1-1));
            clear dhfilter1;

        case 'KF'       %% Kichhoff-Fresnel approximation
            dhfilter = exp(-1i*reconParam.dist(i)*phaseConst*digiHoloTemplKL);
        otherwise
    end
    
    switch( reconParam.outType )
        case 'Intensity'
            recImage=abs(ifft2( dhImgFFT.*dhfilter ));
            idx =  recImage >= 1 ;
            recImage(idx) = 1;
            reconImages(:,:,i) = uint8( recImage(1:holoImgSiz(1), 1:holoImgSiz(2))*255 );
        case 'complex'
            recImage = ifft2( dhImgFFT.*dhfilter );
            reconImages(:,:,i) = recImage(1:holoImgSiz(1), 1:holoImgSiz(2));

    end
end

return;

