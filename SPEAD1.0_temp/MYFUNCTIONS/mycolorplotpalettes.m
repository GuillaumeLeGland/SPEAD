function [bgyrp_dark,bgyrp_soft] = mycolorplotpalettes()

%http://figuredesign.blogspot.com/2012/04/meeting-recap-colors-in-figures.html

bluedark = [33 64 154]/255;
bluesoft = [71 195 211]/255;

greendark = [8 135 67]/255;
greensoft = [173 209 54]/255;

yellow = [255 222 23]/255; 
orange = [241 140 34]/255;

reddark = [127.5 0 0]/255;
redsoft = [237 28 36]/255;

purple = [150 100 155]/255;
pink = [238 132 181]/255;

bgyrp_dark(1,:) = bluedark;
bgyrp_dark(2,:) = greendark;
bgyrp_dark(3,:) = orange;
bgyrp_dark(4,:) = reddark;
bgyrp_dark(5,:) = purple;

bgyrp_soft(1,:) = bluesoft;
bgyrp_soft(2,:) = greensoft;
bgyrp_soft(3,:) = yellow;
bgyrp_soft(4,:) = redsoft;
bgyrp_soft(5,:) = pink;

return
