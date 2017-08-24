function [BARRIER,BARRIER_Coord,area] = ThroatBarrier3D_ADV(v,normal,A,r_cr,areas_v_min, dims, Gpluses, Gminuses)

% -- Computate CANDIDATE throat barrier and area corresponding to a CANDIDATE THROAT -- % 
% 26-connected dilation searching process
% Enforce a "Condition" to find only significant throat in 3D

% Inputs:
% v - coords of "base" voxel center
% normal - local unit normal vector to the plane P(v,t,theta,phi) computed
%          from function tangentDir.m & normalDir.m
% A - the pore-grain binary images (3D structure,255-pore,0-grain)
% r_cr - Applied in "Condition" 
% areas_v_min - "benchmark candidate" throat area in finding Av = min(Ak) to speed up computation

% Outputs: 
% BARRIER - linear indices of candidate barrier; 
% BARRIER_Coord - coords of candidate barrier(consistent with Avizo format)
% area - candidate throat area
% A - pore-grain binary images after marking for barrier voxels (optional output)

% Marking index: marking to speed up searching process (only search and examine unexamined normal void voxel)
% 255 - normal void voxel
% 150 - barrier voxel; 
% 100 - base voxel; 
% 230 - checked non-barrier void voxels

tic;
nx = dims(1);
ny = dims(2);
barrier = sub2ind(dims, v(2)+1, v(1)+1, v(3)+1);
marks = containers.Map('KeyType','uint32','ValueType','logical');
marks(barrier) = true;
BARRIER = barrier;
BARRIER_Coord = v;

% Compute initial throat surface area and effective throat radius
area = polyArea3D_ADV(BARRIER_Coord,v,normal);
r_eff = sqrt(area/pi());

while ~isempty(barrier) && ~isnan(area) && r_eff <= r_cr
    barrier_new = [];
    barrier_coord_new = [];
    for i = 1:length(barrier)
        ind = barrier(i);
        mask_0 = [ind-nx-1, ind-1, ind+nx-1;
                  ind-nx,   ind,   ind+nx;
                  ind-nx+1, ind+1, ind+nx+1];
        mask_n1 = mask_0-nx*ny;
        mask_p1 = mask_0+nx*ny;
        mask_in = cat(3,mask_n1,mask_0,mask_p1);
        mask = maskMod(ind,mask_in, Gpluses, Gminuses);
        for j = 1:size(mask,1)*size(mask,2)*size(mask,3)
            if A(mask(j)) == 255 && ~isKey(marks, mask(j))  % only normal pore voxels (haven't been examined) will be checked
                [I,J,K] = ind2sub(dims, mask(j));
                void_coord = [J-1,I-1,K-1];
                corner_coord = [void_coord(1)-0.5, void_coord(2)-0.5, void_coord(3)-0.5;
                                void_coord(1)-0.5, void_coord(2)+0.5, void_coord(3)-0.5;
                                void_coord(1)+0.5, void_coord(2)-0.5, void_coord(3)-0.5;
                                void_coord(1)+0.5, void_coord(2)+0.5, void_coord(3)-0.5;
                                void_coord(1)-0.5, void_coord(2)-0.5, void_coord(3)+0.5;
                                void_coord(1)-0.5, void_coord(2)+0.5, void_coord(3)+0.5;
                                void_coord(1)+0.5, void_coord(2)-0.5, void_coord(3)+0.5;
                                void_coord(1)+0.5, void_coord(2)+0.5, void_coord(3)+0.5];
                w = zeros(8);
                for k = 1:size(corner_coord,1)
                    vTocorner = corner_coord(k,:) - v; % vector from v to the 8 corners of voxel under examination
                    w(k) = dot(normal,vTocorner); % compute dot products
                end
                % Intersecting test: whether the plane P passes through the voxel under examination
                if (~isempty(find(w>0,1)) && ~isempty(find(w<0,1))) || (~isempty(find(w==0,1)) && length(find(w>=0))==length(w))
                    barrier_new = [barrier_new;mask(j)];
                    barrier_coord_new = [barrier_coord_new; void_coord];
                    %A(mask(j))=150; % mark as barrier voxel
                    marks(mask(j)) = true;
                else
                    %A(mask(j))=230; % mark as non-barrier voxel
                    marks(mask(j)) = true;
                end
            end
        end
    end
    BARRIER = [BARRIER;barrier_new];
    BARRIER_Coord = [BARRIER_Coord;barrier_coord_new];
    area = area + polyArea3D_ADV(barrier_coord_new,v,normal);
    r_eff = sqrt(area/pi());
    % speed up computation by comparing the current computed area with already computed minimal area
    % Apply "Condition" - if r_eff > r_cr, ignore plane P as a possible throat at v
    if ~isnan(areas_v_min) && (area > areas_v_min || r_eff > r_cr)
        BARRIER = nan;
        BARRIER_Coord = nan;
        area = nan;
    end
    barrier = barrier_new;
end

end