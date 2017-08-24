% All image related information and other input data are in this .mat file 
load SEG_RVE_Info.mat

% global variables that have been created contain:
% Gminuses
% Gpluses
% Gxminus
% Gxplus
% Gyminus
% Gyplus
% Gzminus
% Gzplus
% dims
% nx
% ny
% nz
% these variables will be used in different functions 


% SEG_SUB is the one can be splitted into pieces (by row)
tic;
LIST_SUB = ThroatFind3D_ADV(SEG_SUB,Alpha);
t_ThroatFind3D = toc;
toc;
save('ThroatFind3D','-v7.3');