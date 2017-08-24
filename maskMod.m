function mask = maskMod(ind,mask_in)

% Modify masks used in 26-connected dilation searching process
% Consider six boundary planes condition

global Gxplus
global Gyplus
global Gzplus
global Gxminus
global Gyminus
global Gzminus

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