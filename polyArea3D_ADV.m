function area = polyArea3D_ADV(BARRIER_Coord,v,normal)

% -- Compute intersecting area of plane P with barrier voxels -- %
% Complement candidate throat computation in dilation searching
% Each intersecting area is a 3D PLANAR POLYGON 

% Inputs:
% BARRIER_Coord - barrier voxels coords(consistent with Avizo format)
% v - coords of "base" voxel center  
% normal - local unit normal vector to the plane P(v,t,theta,phi) computed
%          from function tangentDir.m & normalDir.m

unit_normal = normal; 


cc = BARRIER_Coord;  
area = 0;
for i = 1:size(cc,1)
    c = cc(i,:); % coords of center of each voxel

    % Determine intesecting points of plane P and each barrier voxel by
    % calculating intesections with all the 12 edges (divided into 3 groups of edges) 
    intersects = [];
    edge_G1 = [c(1)-0.5,c(2)-0.5;
               c(1)+0.5,c(2)-0.5;
               c(1)-0.5,c(2)+0.5;
               c(1)+0.5,c(2)+0.5];
    z = (dot(v,unit_normal)-unit_normal(1)*edge_G1(:,1)-unit_normal(2)*edge_G1(:,2))/unit_normal(3);
    for j = 1:length(z)
        if z(j)>=c(3)-0.5 && z(j)<=c(3)+0.5
            intersects = [intersects;edge_G1(j,:),z(j)];
        end
    end
    edge_G2 = [c(2)-0.5,c(3)-0.5;
               c(2)+0.5,c(3)-0.5;
               c(2)-0.5,c(3)+0.5;
               c(2)+0.5,c(3)+0.5];
    x = (dot(v,unit_normal)-unit_normal(2)*edge_G2(:,1)-unit_normal(3)*edge_G2(:,2))/unit_normal(1);
    for j = 1:length(x)
        if x(j)>=c(1)-0.5 && x(j)<=c(1)+0.5
            intersects = [intersects;x(j),edge_G2(j,:)];
        end
    end
    edge_G3 = [c(1)-0.5,c(3)-0.5;
               c(1)+0.5,c(3)-0.5;
               c(1)-0.5,c(3)+0.5;
               c(1)+0.5,c(3)+0.5];
    y = (dot(v,unit_normal)-unit_normal(1)*edge_G3(:,1)-unit_normal(3)*edge_G3(:,2))/unit_normal(2);
    for j = 1:length(y)
        if y(j)>=c(2)-0.5 && y(j)<=c(2)+0.5
            intersects = [intersects;edge_G3(j,1),y(j),edge_G3(j,2)];
        end
    end

    % Sort all interesection points (vertices of a 3D planar polygon) into counterclockwise orientation 
    % when viewed from the side of plane P pointed to by unit_normal vector
    CENTER = mean(intersects);
    v_intersects = [];
    for j = 1:size(intersects,1)
        v_intersects = [v_intersects;intersects(j,:)-CENTER];
    end
    % "intersects = []" Case: Due to "Tolerance and Stopping Criteria" setting in computing "tangent", 
    % and floating-point accuracy in computing "intersects", this difference usually 
    % corresponds to a very small intersecting area (should have negligible effect on computed area)
    if isempty(intersects)  
        VERT = [];
    elseif ~isempty(intersects)
        VERT = intersects(1,:);
        VERT_pos = [];
        VERT_neg = [];
        v_pos = [];
        v_neg = [];
        for j = 2:size(intersects,1)
            if dot(unit_normal,cross(v_intersects(1,:),v_intersects(j,:)))>=0
                VERT_pos = [VERT_pos;intersects(j,:)];
                v_pos = [v_pos;v_intersects(j,:)];
            else
                VERT_neg = [VERT_neg;intersects(j,:)];
                v_neg = [v_neg;v_intersects(j,:)];
            end
        end
        COSIN = [];
        for j = 1:size(VERT_pos,1)
            COSIN = [COSIN;dot(v_intersects(1,:),v_pos(j,:))/norm(v_intersects(1,:))/norm(v_pos(j,:))];
        end
        [~,I] = sort(COSIN,'descend');
        VERT_pos = VERT_pos(I,:);
        COSIN = [];
        for j = 1:size(VERT_neg,1)
            COSIN = [COSIN;dot(v_intersects(1,:),v_neg(j,:))/norm(v_intersects(1,:))/norm(v_neg(j,:))];
        end
        [~,I] = sort(COSIN,'ascend');
        VERT_neg = VERT_neg(I,:);

        VERT = [VERT;VERT_pos;VERT_neg]; % 3D polygon vertices information
    end

    % Compute 3D polygon area 
    % (ref:http://geomalgorithms.com/a01-_area.html#3D%20Polygons, 
    % http://stackoverflow.com/questions/12642256/python-find-area-of-polygon-from-xyz-coordinates)
    if size(VERT,1)<3
        AREA = 0;
    else
        sumprod = [0,0,0];
        for j = 1:size(VERT,1)
            v1 = VERT(j,:);
            if j == size(VERT,1)
                v2 = VERT(1,:);
            else
                v2 = VERT(j+1,:);
            end
            crprod = cross(v1,v2);
            sumprod = sumprod + crprod;
        end
        AREA = dot(unit_normal,sumprod)/2;
    end
    area = area + AREA;
end

end