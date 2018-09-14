%%%------------------------------------
% node ID to coordinate
%%%------------------------------------

function [coord] = IDtoCoord(nodeID,nx,ny,dimx,dimy,dimz)

z = (floor((nodeID-1)/(nx*ny))+1);
y = 1+floor((nodeID-.01-(z-1)*nx*ny)/nx);
x = nodeID-((z-1)*nx*ny)-((y-1)*nx);
z=z*dimz; y=y*dimy; x=x*dimx;

coord=[x,y,z];