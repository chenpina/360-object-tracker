function [result] = get_sceen(oriImage, cr, cc, rh, rw, scale, draw)
% [result] = getSceen(origin_image, center_row, center_col, result_hight, result_width, scale_factor, draw_flag)
% size(result) = [rh/scale(1),rw/scale(2),size(oriImage,3)]
%% calculate parameters(image size, sphere radius)
if isa(oriImage,'uint8'), oriImage = double(oriImage)/255; end
[ih,iw,id] = size(oriImage);
% max_edge = ih*sqrt(2)/pi;
% if rh > max_edge, rh = floor(max_edge); end
% if rw > max_edge, rw = floor(max_edge); end
r = ih/pi;
if nargin < 6, scale = [1,1]; end
if length(scale) == 1, scale = [scale,scale]; end
rh = round(rh/scale(1));
rw = round(rw/scale(2));

%% make center point to image center colomn
padImage = padarray(oriImage, [0 iw], 'circular', 'post');
w_c = floor(cc+iw/2+1);
if w_c > iw, w_c = w_c - iw; end
image = padImage(:,w_c:w_c+iw-1,:);
clear padImage

%% calculate degree of eye point to result sceen's points
[c_vec(1,1,1),c_vec(1,1,2),c_vec(1,1,3)] = sph2cart(pi,pi/2-cr*pi/ih,r);
lr = cross(c_vec,reshape([0,0,1],[1,1,3]));
lr = lr / norm(lr(:));
ub = cross(c_vec,lr);
ub = ub / norm(ub(:));
p = repmat(c_vec,[rh,rw]) + repmat((1-floor(rh*scale(1)/2):scale(1):ceil(rh*scale(1)/2))',[1,rw,3]).*repmat(ub,[rh,rw])...
                             + repmat(1-floor(rw*scale(2)/2):scale(2):ceil(rw*scale(2)/2),[rh,1,3]).*repmat(lr,[rh,rw]);
[theta,phi,~] = cart2sph(p(:,:,1),p(:,:,2),p(:,:,3));
theta = mod(theta, 2*pi);

%% change degree to position on image
pos_r = (0.5 - phi/pi) * ih;
pos_r(pos_r<1) = 1;
pos_c = (2*pi - theta)/(2*pi) * iw + 1;

%% interpolation and get new color for each points on result sceen
% result = zeros(rh,rw,id);
% for i = 1:rh
%     for j = 1:rw
%         for k = 1:id
%             result(i,j,k) = interpn(image,pos_r(i,j),pos_c(i,j),k);
%         end
%     end
% end
interp_r = cat(4,floor(pos_r),floor(pos_r)+1,floor(pos_r),floor(pos_r)+1); % rh x rw x 1 x 4
interp_c = cat(4,floor(pos_c),floor(pos_c),floor(pos_c)+1,floor(pos_c)+1); % rh x rw x 1 x 4
interp_r(interp_r>ih) = ih;
interp_c(interp_c>iw) = iw;
interp_r = repmat(interp_r,[1,1,id]); % rh x rw x id x 4
interp_c = repmat(interp_c,[1,1,id]); % rh x rw x id x 4
idx_d = ones(rh,rw,id,4);
for i = 2:id, idx_d(:,:,i,:) = i*idx_d(:,:,i,:); end
color = reshape(image(interp_r(:)+(interp_c(:)-1)*ih+(idx_d(:)-1)*ih*iw), [rh,rw,id,4]);
diff_r = pos_r - interp_r(:,:,1,1); % rh x rw
diff_c = pos_c - interp_c(:,:,1,1); % rh x rw
ratio(:,:,1,1) = (1-diff_r) .* (1-diff_c);
ratio(:,:,1,2) = diff_r .* (1-diff_c);
ratio(:,:,1,3) = (1-diff_r) .* diff_c;
ratio(:,:,1,4) = diff_r .* diff_c;
ratio = repmat(ratio,[1,1,id]); % rh x rw x id x 4
result = sum(color.*ratio,4); % rh x rw x id = sum(rh x rw x id x 4, 4);

%% show result
if nargin > 6
    imshow(result)
end
end