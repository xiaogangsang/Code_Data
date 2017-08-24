function normvec = normalDir(tangent)

% -- Define a sequence of unit normal vectors for the plane P including theta and phi -- %
% all angles computed are in degrees
% 
% INPUT:
% tangent - unit vector for local tangent direction to the skeleton computed from function - tangentDir.m
% 
% OUTPUT:
% normvec - [theta, phi, normvec x, normvec y, normvec z]

tx = tangent(1);
ty = tangent(2);
tz = tangent(3);

% thetat range - [0,180]; phit range - [0,360)
thetat = acosd(tz/norm(tangent));
phit = atand(ty/tx); % "atand" return phit in [-90,90]
if tx > 0 && ty >= 0  %[0,90)
    phit = phit;
elseif tx < 0 && ty >= 0 %(90,180] 
    phit = phit + 180; 
elseif tx < 0 && ty < 0 %(180,270)
    phit = phit + 180;
elseif tx > 0 && ty < 0 %(270,360)
    phit = phit + 360;
elseif tx == 0 && ty > 0 %90
    phit = phit;
elseif tx == 0 && ty < 0 %270
    phit = phit + 360;
elseif tx == 0 && ty == 0 % given NAN
    phit = 0; % thetat must be 0
end

% form rotation matrix to rotate "normvec" along with "tangent"
Ry = [cosd(thetat),0,sind(thetat);
      0,1,0;
      -sind(thetat),0,cosd(thetat)];
Rz = [cosd(phit),-sind(phit),0;
      sind(phit),cosd(phit),0;
      0,0,1];

% assume "tangent" as the zenith axis in a spherical coordinate system
% centered at the base point, "norm_t" is specified by theta and phi wrt. "tangent", 
% "normvec" is computed by premultipling "norm_t" by Ry and Rz to get
% the unit vector wrt. XYZ coords
% normvec = [theta, phi, normvec x, normvec y, normvec z]

%%%%%%%%%%%%SUBJECT TO CHANGE%%%%%%%%%%%%
theta = 5:5:45;
phi = 0:45:315;

% theta = 15:15:45;
% phi = 0:45:315;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% one particular case is theta = 0 & phi = 0 --> n = t
normvec = zeros(length(theta)*length(phi)+1,5);
angles = [];
for i = 1:length(theta)
    for j = 1:length(phi)
        angles = [angles;theta(i),phi(j)];
    end
end
normvec(2:size(normvec,1),1:2) = angles;

for i = 1:size(normvec,1)
    theta = normvec(i,1);
    phi = normvec(i,2);
    norm_t = [sind(theta)*cosd(phi);
              sind(theta)*sind(phi);
              cosd(theta)];  % column vector
    norm_t_rot = Rz*Ry*norm_t;
    normvec(i,3) = norm_t_rot(1);
    normvec(i,4) = norm_t_rot(2);
    normvec(i,5) = norm_t_rot(3);
end
    
end