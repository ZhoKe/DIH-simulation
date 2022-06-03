function [Xp, Yp, Zp] = P_seeding(xlim, ylim, zlim, Np)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here


Lx = xlim(2) - xlim(1);
Ly = ylim(2) - ylim(1);
Lz = zlim(2) - zlim(1);

Xp = xlim(1) + (Lx * rand([Np,1]));
Yp = ylim(1) + (Ly * rand([Np,1]));
Zp = zlim(1) + (Lz * rand([Np,1]));

end

