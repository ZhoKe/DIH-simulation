function [digiHoloTemplKL] = createHoloTempl( imgSiz )
% function for creating filter
%   imgSiz - [col, row]
colNum = 2^ceil( log(imgSiz(1))/log(2) );   % number of column
rowNum = 2^ceil( log(imgSiz(2))/log(2) );   % number of row

templK = zeros( [rowNum, colNum] );
templL = zeros( [rowNum, colNum] );

for i=1:colNum    % column
    templL(:,i) = i;
end
templL = fftshift( templL- (colNum/2+1), 2 );
templL = templL/colNum;

for i=1:rowNum
    templK(i,:) = i;
end
templK = fftshift( templK- (rowNum/2+1), 1 );
templK = templK/rowNum;

digiHoloTemplKL = templK.^2+templL.^2;
return;
