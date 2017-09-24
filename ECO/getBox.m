function [result] = getBox(oriImage, cr, cc, rh, rw, draw)
% [result] = getBox(origin_image, center_row, center_col, result_hight, result_width, draw_flag)
% result is 2x2*(rh+rw) matrix
%% calculate parameters(image size, sphere radius)
[ih,iw,id] = size(oriImage);
% max_edge = ih*sqrt(2)/pi;
% if rh > max_edge, rh = floor(max_edge); end
% if rw > max_edge, rw = floor(max_edge); end
r = ih/pi;

%% calculate degree of eye point to result box edge's points
[c_vec(1,1),c_vec(2,1),c_vec(3,1)] = sph2cart(pi,pi/2-cr*pi/ih,r);
lr = cross(c_vec,[0;0;1]);
lr = lr / norm(lr(:));
ub = cross(c_vec,lr);
ub = ub / norm(ub(:));
p = repmat(c_vec,[1,2*(rh+rw)]) + repmat([1-rh/2:rh/2,(rh/2)*ones(1,rw),rh/2:-1:1-rh/2,(1-rh/2)*ones(1,rw)],[3,1]).*repmat(ub,[1,2*(rh+rw)])...
                              + repmat([(1-rw/2)*ones(1,rh),1-rw/2:rw/2,(rw/2)*ones(1,rh),rw/2:-1:1-rw/2],[3,1]).*repmat(lr,[1,2*(rh+rw)]);
[theta,phi,~] = cart2sph(p(1,:),p(2,:),p(3,:));
theta = mod(theta, 2*pi);

%% change degree to position on image
pos_r = (0.5 - phi/pi) * ih;
pos_r(pos_r<1) = 1;
pos_c = (2*pi - theta)/(2*pi) * iw + 1;
pos_c = mod(pos_c + iw/2 + cc, iw);
pos_c(pos_c==0) = iw;
result = [pos_r;pos_c];

%% show result
if nargin > 5
    imshow(oriImage)
    hold on
    plot(result(2,:),result(1,:),'r.')
    plot(cc,cr,'g.')
    hold off
end
end