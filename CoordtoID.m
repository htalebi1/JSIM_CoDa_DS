%%%------------------------------------
% coordinate to node ID
%%%------------------------------------

function [nodeID,out_win] = CoordtoID(coord,nx,ny,nz,dimx,dimy,dimz)

coord=coord./repmat([dimx,dimy,dimz],size(coord,1),1);
nodeID = (coord(:,3)-1)*nx*ny + (coord(:,2)-1)*nx + coord(:,1);
nodeID(nodeID<1)=1;
nodeID(nodeID>nx*ny*nz)=1;
out_win = [find(coord(:,1)<0.1)',find(coord(:,2)<0.1)',find(coord(:,3)<0.1)',...
           find(coord(:,1)>nx)',find(coord(:,2)>ny)',find(coord(:,3)>nz)',find(nodeID<1)',find(nodeID>nx*ny*nz)'];

