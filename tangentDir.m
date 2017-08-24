function tangent = tangentDir(v,p)

% Compute local tangent direction to the skeleton based on each "base" voxel
% Use a least-squares fitting of base voxel and its neighboring six skeleton voxels
% v is a vector containing the coordinates of center of base voxel
% p is a matrix containing the coordinates of neghboting voxels' centers and v

syms a b c

t = [a,b,c]; % vector form of the tangent direction 

di = [];
for i = 1:size(p,1)
    vvi = p(i,:)-v;
    di_square = sum(cross(vvi,t).*cross(vvi,t))./sum(t.*t);
    di = [di;sqrt(di_square)];
end

guess = p(size(p,1),:)-p(1,:); % set starting point 
%options = optimoptions('lsqnonlin','TolFun',1e-10,'TolX',1e-10,'Display','off'); 
options = optimoptions('lsqnonlin','Display','off');
fun = @(x)double(subs(di,[a,b,c],x));
x = lsqnonlin(fun,guess,[],[],options);

tangent = x;
% normalize to get unit local tangent vector
tangent = tangent/norm(tangent);


end
