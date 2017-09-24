function [img_files, pos, target_sz, ground_truth, video_path] = load_video_info_360(video_path)

ground_truth = dlmread([video_path '/init.txt']);

frames = size(ground_truth, 1);
seq.init_rect = ground_truth(1,:);

%set initial position and size
target_sz = [ground_truth(1,3), ground_truth(1,4)];
pos = [ground_truth(1,1), ground_truth(1,2)];

img_path = [video_path '/img/'];

%see if they are in the 'img' subfolder or not
if exist([img_path num2str(1, '%04i.png')], 'file')
    img_files = num2str((1:frames)', [img_path '%04i.png']);
elseif exist([img_path num2str(1, '%04i.jpg')], 'file')
    img_files = num2str((1:frames)', [img_path '%04i.jpg']);
elseif exist([img_path num2str(1, '%04i.bmp')], 'file')
    img_files = num2str((1:frames)', [img_path '%04i.bmp']);
else
    error('No image files to load.')
end

%list the files
img_files = cellstr(img_files);

end

