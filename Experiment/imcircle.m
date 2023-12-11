%__________________________________________________________________________
function [mask,ind] = imcircle(ix,iy,cx,cy,r)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-Define a total field (circular mask) in which dots could be placed.

%-mask: a matrix in which all cell inside the circle is one
%-ind: index of the cell equals to one
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[x,y]=meshgrid(-(cx-1):(ix-cx),-(cy-1):(iy-cy)); 
mask=((x.^2+y.^2)<=r^2);
ind = find(mask~=0);

end
