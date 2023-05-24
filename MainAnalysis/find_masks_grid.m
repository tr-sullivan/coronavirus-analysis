function [masks_temp, masks_bg] = find_masks_grid(rgb1, int_corr, radius_range, rnd, array_size, shifts)
% finds masks using a semi-interactive procedure
% parameter inputs:
% rgb1 is the image to segment
% int_corr is used to rescale the input image so the range 0-int_corr is
% used
% radius_range limits the range of possible radii extracted by the circle
% detection algorithm
% rnd is used to round so spots within the same row/column are grouped
% shifts specifies the x shift along the grid between each row

nRow = array_size(1);
nCol = array_size(2);

figure();imshow(rgb1);
K = imadjust(rgb1,[0 int_corr],[]);

% restrict circle detection to only within a certain range of the image
[hei, wid] = size(K);
vertthresh = [0.00 0.00];
horizthresh = [0.01 0.01];

hthresh = [floor(vertthresh(1)*hei), hei - floor(vertthresh(2)*hei)];
wthresh = [floor(horizthresh(1)*wid), wid - floor(horizthresh(2)*wid)];

% initial sensitivity and threshold for circle detection
se(1) = 0.93; se(2) = 0.01;
[centers, radii] = imfindcircles(K,radius_range,'ObjectPolarity','bright', 'Sensitivity',se(1), 'EdgeThreshold', se(2));
imshow(K)
line([wthresh(1) wthresh(1)], [hthresh(1), hthresh(2)]);
line([wthresh(2) wthresh(2)], [hthresh(1), hthresh(2)]);
line([wthresh(1) wthresh(2)], [hthresh(1), hthresh(1)]);
line([wthresh(1) wthresh(2)], [hthresh(2), hthresh(2)]);

% restrict to within specified region and add background ROI
centers = centers(centers(:,1) > wthresh(1) & centers(:,1) < wthresh(2) & centers(:,2) > hthresh(1) & centers(:,2) < hthresh(2),:);
radii = radii(centers(:,1) > wthresh(1) & centers(:,1) < wthresh(2) & centers(:,2) > hthresh(1) & centers(:,2) < hthresh(2),:);

% visualize circles, add in radii as separate column and subtract so
% first two columns represent corner

h = viscircles(centers,radii, 'EdgeColor','w', 'LineWidth', 5);
centers(:,3) = radii(:,1);

corners = zeros(4,2);
fprintf('Click on center of top-left corner spot.\n');
corners(1,:) = ginput(1);
fprintf('Click on center of top-right corner spot.\n');
corners(2,:) = ginput(1);
fprintf('Click on center of bottom-left corner spot.\n');
corners(3,:) = ginput(1);
fprintf('Click on center of bottom-right corner spot.\n');
corners(4,:) = ginput(1);

%estimate grid spacing
Dx1 = corners(2,1) - corners(1,1);
Dx2 = corners(4,1) - corners(3,1);
Dy1 = corners(3,2) - corners(1,2);
Dy2 = corners(4,2) - corners(2,2);

Dx = 0.5*(Dx1 + Dx2);
Dy = 0.5*(Dy1 + Dy2);

Dtheta1 = corners(2,2) - corners(1,2);
Dtheta2 = corners(4,2) - corners(3,2);
Dtheta = 0.5*(Dtheta1 + Dtheta2);

dx = Dx/(array_size(2) - 1);
dy = Dy/(array_size(1) - 1);
dtheta = Dtheta/(array_size(2) - 1);

%find center point of each grid element
grid = zeros(array_size(1), array_size(2), 2);
grid(1,1,:) = corners(1,:);
for i = 1:array_size(1)
    x = corners(1,1) + dx*sum(shifts(1:(i-1)));
    y = corners(1,2) + dy*(i-1);
    for j = 1:array_size(2)
        grid(i,j,1) = x + (j-1)*dx - dtheta*(i-1);
        grid(i,j,2) = y + dtheta*(j-1);
    end
end

selected_spots = zeros(array_size(1), array_size(2), 3);
rnz = [];
xshift_real = [];
yshift_real = [];
for i = 1:array_size(1)
    for j = 1:array_size(2)
        point = grid(i,j,:);
        candidates_x = abs(centers(:,1) - point(1)) < .5*dx;
        candidates_y = abs(centers(:,2) - point(2)) < .5*dy;
        if ~isempty(rnz)
            candidates_r = abs(centers(:,3) - median(rnz)) < 5;
            candidates = candidates_x & candidates_y & candidates_r;
        else
            candidates = candidates_x & candidates_y;
        end
        if ~any(candidates)
            selected_spots(i,j,1) = point(1);
            selected_spots(i,j,2) = point(2);
        elseif sum(sum(candidates)) == 1
            id = find(candidates);
            selected_spots(i,j,1) = centers(id,1);
            selected_spots(i,j,2) = centers(id,2);
            selected_spots(i,j,3) = centers(id,3);
            rnz = [centers(id,3), rnz];
            xshift_real = [xshift_real, point(1) - centers(id,1)];
            yshift_real = [yshift_real, point(2) - centers(id,2)];
        else
            candidate_coords = find(candidates);
            dist_x = abs(centers(candidates,1) - point(1));
            dist_y = abs(centers(candidates,2) - point(2));
            distance = sqrt(dist_x .^ 2 + dist_y .^ 2);
            id = candidate_coords(find(distance == min(distance)));
            selected_spots(i,j,1) = centers(id,1);
            selected_spots(i,j,2) = centers(id,2);
            selected_spots(i,j,3) = centers(id,3);
            rnz = [centers(id,3),rnz];
            xshift_real = [xshift_real, point(1) - centers(id,1)];
            yshift_real = [yshift_real, point(2) - centers(id,2)];
        end
        
    end
end

%if any radii are zero (i.e. no corresponding found circles)
medrad = median(rnz);
for i = 1:array_size(1)
    for j = 1:array_size(2)
        if selected_spots(i,j,3) == 0
            selected_spots(i,j,3) = medrad;
        end
    end
end
disc_radius = min(dx, dy)/2 - medrad;

%plotting spots
close all;
figure();
imshow(K)
hold on

% final check of the spots and making masks
cen = selected_spots(:,:,1:2);
cen(:,:,1) = cen(:,:,1) - selected_spots(:,:,3);
cen(:,:,2) = cen(:,:,2) - selected_spots(:,:,3);
cen(:,:,3) = selected_spots(:,:,3);

for i = 1:array_size(1)
    for j = 1:array_size(2)
        % plot detected ROI
        %text(cen(i,j,1), cen(i,j,2) ,num2str((i-1)*array_size(2)+j), 'color', 'w', 'fontsize', 1); 
        h = imellipse(gca, [cen(i,j,1),cen(i,j,2), cen(i,j,3)*2, cen(i,j,3)*2]);
        setColor(h,'g');
        % plot background ROI
        h2 = imellipse(gca, [cen(i,j,1)-disc_radius,cen(i,j,2)-disc_radius,...
            (cen(i,j,3)+disc_radius)*2, (cen(i,j,3)+disc_radius)*2]);
        setColor(h2,'y');
        %position = wait(timer);
        text(cen(i,j,1)-20, cen(i,j,2)-5 ,num2str((i-1)*array_size(2)+j), 'color', 'w', 'fontsize', 12); 
        drawnow; 
        
        % make mask
        masks_temp{(i-1)*array_size(2)+j,1} = createMask(h);
        masks_bg{(i-1)*array_size(2)+j,1} = analyze_background(reshape(cen(i,j,1:3),[1,3]) + [cen(i,j,3),cen(i,j,3),0], size(rgb1), disc_radius);
    end
end
saveas(gcf, ['Chip-1-mask'] , 'pdf');
% end checking

end