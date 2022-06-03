function [Up_Holo] = Contrast_recov ( HighC_Holo, LowC_Holo)

HighC_Holo = double(HighC_Holo);
LowC_Holo = double( LowC_Holo);

I_H1 = mean( HighC_Holo(:));
I_H2 = mean( LowC_Holo(:));
S_H1 = std( HighC_Holo(:));
S_H2 = std( LowC_Holo(:));
Up_Holo = I_H1 + (LowC_Holo - I_H2).*(S_H1/S_H2);

end