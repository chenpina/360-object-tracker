
% This demo script runs the ECO tracker with deep features on the
% included "Crossing" video.

% Add paths
setup_paths();

% Load video information
datafilename = 'aZQ1OBpTCvY';
fprintf('%s\n', datafilename);
video_path = ['sequences/',datafilename];
[seq, ground_truth] = load_video_info(video_path);

% Run ECO
results = testing_ECO_360(seq);

% show result for 360
s_pathname = ['sequences/',datafilename,'/'];
im = imread(num2str(1,[s_pathname,'img/%04i.png']));
scale = 1;
while size(im,2)/scale > 1600, scale = scale*2; end
imshow(im);
set(gcf, 'pos', [100,100,size(im,2)/scale,size(im,1)/scale]);
if size(im,1)/scale==720
    set(gcf, 'pos', [100,100,size(im,2)/scale*3/4,size(im,1)/scale*3/4]);
end
for i=1:size(results.res,1)
    res = results.res(i,:);
    im = imread(num2str(i,[s_pathname,'img/%04i.png']));
    getBox(im, res(1), res(2), res(3), res(4), []);
    set(gca, 'pos', [0,0,1,1]);
end
