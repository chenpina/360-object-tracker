function [result] = p2Dto360(oriImage, cr, cc, rh, rw, pr, pc, sh, sw)
% [result] = getBox(origin_image, center_row, center_col, result_hight, result_width, point_row, point_col, scale_factor)
%% calculate parameters(image size, sphere radius)
[ih,iw,~] = size(oriImage);
% max_edge = ih*sqrt(2)/pi;
% if rh > max_edge, rh = floor(max_edge); end
% if rw > max_edge, rw = floor(max_edge); end
r = ih/pi;
if nargin < 8
    sh = rh; sw = rw;
elseif nargin < 9
    sw = sh;
end
scale = floor([sh/rh,sw/rw]);
scale(scale<1) = 1;
sh = 2*round(sh/scale(1)/2);
sw = 2*round(sw/scale(2)/2);

%% calculate degree of eye point to result box edge's points
[c_vec(1,1),c_vec(2,1),c_vec(3,1)] = sph2cart(pi,pi/2-cr*pi/ih,r);
lr = cross(c_vec,[0;0;1]);
lr = lr / norm(lr(:));
ub = cross(c_vec,lr);
ub = ub / norm(ub(:));
pr = (pr-1)*sh/rh + 1;
pc = (pc-1)*sw/rw + 1;
p = c_vec + (1+(pr-sh/2-1)*scale(1))*ub + (1+(pc-sw/2-1)*scale(2))*lr;
[theta,phi,~] = cart2sph(p(1),p(2),p(3));
theta = mod(theta, 2*pi);

%% change degree to position on image
pos_r = (0.5 - phi/pi) * ih;
pos_r(pos_r<1) = 1;
pos_c = (2*pi - theta)/(2*pi) * iw + 1;
pos_c = mod(pos_c + iw/2 + cc, iw);
pos_c(pos_c==0) = iw;
result = [pos_r,pos_c];
end