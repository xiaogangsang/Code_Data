function LIST = ThroatFind3D_ADV(Seg,Alpha)

% -- Independently compute the CROSS-SECTION LISTS for Significant Pore Paths -- %

% INPUTS:
% Seg: All/Sub-Set of SEG (SEG: SPP for 3D whole structure)
% Alpha : 'Pore.mat', pore-grain images (3D structure,255-pore,0-grain)
%
% OUTPUTS:
% LIST = {(1)-Seg ID, (2)-Nodal_1, (3)-Nodal_2, (4)-{LIST_seg}, (5)-Current throat area}
% LIST_seg = {(1)-Base voxel ID (Point_ID), 
%             (2)-angle theta for Candidate throat, 
%             (3)-angle phi for Candidate throat,
%             (4)-Candidate throat area Av
%             (5)-Candidate barrier linear indices,
%             (6)-Candidate barrier Coords(X,Y,Z)}
%
% LIST_seg is sorted based on "Candidate throat area", Current throat area = min(Candidate throat area)
% LIST is sorted based on "Current throat area"

LIST =cell(size(Seg,1),5);
for i = 1:size(Seg,1)
    i
    Seg_index = i;
    POINTinfo = Seg{Seg_index,6};
    
    % Do not consider points corresponding to two end nodes of segments (conceptually nodal pore centers)
    LIST_seg = cell(size(POINTinfo,1)-2,6);
    for j = 2:size(POINTinfo,1)-1
        v_index = j;
        v = cell2mat(POINTinfo(v_index,2:4));  % "base" voxel coords
        % Construct vector p containing v and its six neighboring skeleton voxels
        % Involve seven skeleton voxels
        if size(POINTinfo,1) <= 7
            p = cell2mat(POINTinfo(:,2:4));
        elseif size(POINTinfo,1) > 7
            % Consider 4 special cases where the six neighboring voxels are not symmetric 
            if j==2
                p_index = v_index-1:v_index+5;
            elseif j==3
                p_index = v_index-2:v_index+4;
            elseif j==size(POINTinfo,1)-2
                p_index = v_index-4:v_index+2;
            elseif j==size(POINTinfo,1)-1
                p_index = v_index-5:v_index+1;
            else
                p_index = v_index-3:v_index+3; 
            end
            p = cell2mat(POINTinfo(p_index,2:4));
        end
        % Compute local tangent direction t to the skeleton
        tangent = tangentDir(v,p);
        % Compute local unit normal vector n to the plane P(v,t,theta,phi)
        normvec = normalDir(tangent);
        % Extract r_cr = min(Ls,Le) from POINTinfo 
        r_cr = POINTinfo{j,7};
        % Search for a minimal area throat /barrier candidate at v using "normvec"
        % Av = min(Ak), Ak is for different theta and phi
        areas_v = cell(size(normvec,1),1);
        barriers_v = cell(size(normvec,1),1);
        barriers_coord_v = cell(size(normvec,1),1);
        areas_v_min = nan; % speed up computation: compute Av as iterating k
        for k = 1:size(normvec,1)
            normal = normvec(k,3:5);
            [BARRIER_v,BARRIER_Coord_v,area_v] = ThroatBarrier3D_ADV(v,normal,Alpha,r_cr,areas_v_min);
            areas_v{k} = area_v;
            barriers_v{k} = BARRIER_v;
            barriers_coord_v{k} = BARRIER_Coord_v;
            
            areas_v_min = min(cell2mat(areas_v));
        end
        LIST_seg_v = [num2cell(normvec(:,1:2)),areas_v,barriers_v,barriers_coord_v];
        [~,I] = sort(cell2mat(LIST_seg_v(:,3)),'ascend');
        LIST_seg_v = LIST_seg_v(I,:);
        % Create a cross-section list for each SPP in Seg
        LIST_seg{j-1,1} = POINTinfo{j,1};
        LIST_seg(j-1,2:6) = LIST_seg_v(1,:); % info corresponding to Av 
    end
    % Sort LIST_seg based on "Candidate throat area"
    [~,I] = sort(cell2mat(LIST_seg(:,4)),'ascend');
    LIST_seg = LIST_seg(I,:);
    % Create LIST for all SPPs in Seg
    LIST{i,1}=Seg{i,1};
    LIST{i,2}=Seg{i,3};
    LIST{i,3}=Seg{i,4};
    LIST{i,4}=LIST_seg;
    LIST{i,5}=LIST_seg{1,4};
end
% Sort LIST based on "Current throat area"
[~,I] = sort(cell2mat(LIST(:,5)),'ascend');
LIST = LIST(I,:);

end