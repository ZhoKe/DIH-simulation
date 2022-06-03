function [Lz1,Lz,Lz2] = MT_DOFestimator_Fn(HoloSynthLz,rec_imgSynth,ImgSiz, Dp, Lambda, Reso, Z_depth, LimPer);


% DOF or Lz estimator:


frames = 1:length(Z_depth);
Circ_R = Dp/2;
xx = ImgSiz(1)/2;
yy = ImgSiz(2)/2;
Cent_I = 0;
Bckg_I = 120;
[Ap] = Circ_Aperture_Particle(ImgSiz(1),ImgSiz(2), Circ_R, Cent_I, Bckg_I,xx,yy);
[Holo] = MT_HoloGen_Fn2(Lambda,Reso,Ap,mean(Z_depth),0);
imwrite(uint8(abs(Holo)),sprintf( HoloSynthLz), 'tif', 'Compression', 'none' );
[rec] = MT_HoloRecon_Fn( 1,Lambda,Reso,HoloSynthLz,Z_depth, rec_imgSynth );
figure,
del = 60;
sig = zeros(2*del+1,1);
s=1;
for iz = round(mean(frames))-del:round(mean(frames))+del
    zz = iz
    a = imcomplement(imread(sprintf( rec_imgSynth,1,iz)));
    sig(s,1) = double(a(yy,xx));
    s = s+1;
end

ss = (sig-min(sig))./(max(sig)-min(sig));
plot(ss,'r')
IDs = find (ss > LimPer);
Lz1 = max(IDs)-min(IDs)
Lz2 = round(((Dp*Reso)^2)/Lambda)
Lz = min(Lz1,Lz2)
end



