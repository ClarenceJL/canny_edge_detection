%% test
% run cannyEdge for all images

imgpath = '../TestScript_P1/train_images_P1/*.jpg';
images = dir(imgpath);

for i = 1:length(images)
    I = imread(['../TestScript_P1/train_images_P1/' images(i).name]);
    E = cannyEdge(I);
    print(h,'-djpeg',[images(i).name(1:end-4) '_edge.jpg']);
    close all
    %save([['../TestScript_P1/train_images_P1/' images(i).name(1:end-4) '_edge.mat']],'E');
end