function mask = maskMod(ind, mask_in, Gpluses, Gminuses)

% Modify masks used in 26-connected dilation searching process
% Consider six boundary planes condition

Gxplus = Gpluses{1, 1};
Gyplus = Gpluses{1, 2};
Gzplus = Gpluses{1, 3};
Gxminus = Gminuses{1, 1};
Gyminus = Gminuses{1, 2};
Gzminus = Gminuses{1, 3};

if ~isempty(find(ind==Gxplus,1))
    mask_in(3,:,:)=[];
end
if ~isempty(find(ind==Gxminus,1))
    mask_in(1,:,:)=[];
end
if ~isempty(find(ind==Gyplus,1))
    mask_in(:,3,:)=[];
end
if ~isempty(find(ind==Gyminus,1))
    mask_in(:,1,:)=[];
end
if ~isempty(find(ind==Gzplus,1))
    mask_in(:,:,3)=[];
end
if ~isempty(find(ind==Gzminus,1))
    mask_in(:,:,1)=[];
end

mask = mask_in;

end