function [BARRIER,BARRIER_Coord,area] = ThroatBarrier3D_ADV(v,normal,A,r_cr,areas_v_min)

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

global nx
global ny
global nz
global X
global Y
global Z

% Marking index: marking to speed up searching process (only search and examine unexamined normal void voxel)
% 255 - normal void voxel
% 150 - barrier voxel; 
% 100 - base voxel; 
% 230 - checked non-barrier void voxels

barrier = find(X==v(1)&Y==v(2)&Z==v(3));
A(barrier)=150;
BARRIER = barrier;
BARRIER_Coord = v;

% Compute initial throat surface area and effective throat radius
area = polyArea3D_ADV(BARRIER_Coord,v,normal);
r_eff = sqrt(area/pi());


while ~isempty(barrier) && ~isnan(area) && r_eff <= r_cr
    barrier_new = [];
    BARRIER_Coord_new = [];
    area_new = [];
    for i = 1:length(barrier)
        ind = barrier(i);
        mask_0 = [ind-nx-1, ind-1, ind+nx-1;
                  ind-nx,   ind,   ind+nx;
                  ind-nx+1, ind+1, ind+nx+1];
        mask_n1 = mask_0-nx*ny;
        mask_p1 = mask_0+nx*ny;
        mask_in = cat(3,mask_n1,mask_0,mask_p1);
        mask = maskMod(ind,mask_in);
        for j = 1:size(mask,1)*size(mask,2)*size(mask,3)
            if A(mask(j)) == 255  % only normal pore voxels (haven't been examined) will be checked
                void_coord = [X(mask(j)),Y(mask(j)),Z(mask(j))];
                corner_coord = [void_coord(1)-0.5, void_coord(2)-0.5, void_coord(3)-0.5;
                                void_coord(1)-0.5, void_coord(2)+0.5, void_coord(3)-0.5;
                                void_coord(1)+0.5, void_coord(2)-0.5, void_coord(3)-0.5;
                                void_coord(1)+0.5, void_coord(2)+0.5, void_coord(3)-0.5;
                                void_coord(1)-0.5, void_coord(2)-0.5, void_coord(3)+0.5;
                                void_coord(1)-0.5, void_coord(2)+0.5, void_coord(3)+0.5;
                                void_coord(1)+0.5, void_coord(2)-0.5, void_coord(3)+0.5;
                                void_coord(1)+0.5, void_coord(2)+0.5, void_coord(3)+0.5];
                w = [];
                for k = 1:size(corner_coord,1)
                    vTocorner = corner_coord(k,:) - v; % vector from v to the 8 corners of voxel under examination
                    w = [w;dot(normal,vTocorner)]; % compute dot products
                end
                % Intersecting test: whether the plane P passes through the voxel under examination
                if (~isempty(find(w>0,1)) && ~isempty(find(w<0,1))) || (~isempty(find(w==0,1)) && length(find(w>=0))==length(w))
                    barrier_new = [barrier_new;mask(j)];
                    BARRIER_Coord_new = [BARRIER_Coord_new; void_coord];
                    A(mask(j))=150; % mark as barrier voxel
                else
                    A(mask(j))=230; % mark as non-barrier voxel
                end
            end
        end
    end
    BARRIER = [BARRIER;barrier_new];
    BARRIER_Coord = [BARRIER_Coord;BARRIER_Coord_new];
    area_new = polyArea3D_ADV(BARRIER_Coord_new,v,normal);
    area = area + area_new;
    % speed up computation by comparing the current computed area with already computed minimal area
    if ~isnan(areas_v_min) && area > areas_v_min
        BARRIER = nan;
        BARRIER_Coord = nan;
        area = nan;
    else
        % Compute effective throat radius
        r_eff = sqrt(area/pi());
        % Apply "Condition" - if r_eff > r_cr, ignore plane P as a possible throat at v
        if r_eff > r_cr
            BARRIER = nan;
            BARRIER_Coord = nan;
            area = nan;
        end
    end
    barrier = barrier_new;
end

% mark "base" voxel for distinguish
if ~isnan(BARRIER)
    A(BARRIER(1))=100; 
end


end