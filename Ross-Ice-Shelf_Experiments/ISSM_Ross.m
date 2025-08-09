
cd '/Users/rishi/Desktop/ISSM_Project_Ross'

% Load the shapefile
shape = shaperead('Ross_Ice_shelf_Antarctica.shp');

% Extract the coordinates of the polygon
x = shape.X;  % Longitude coordinates
y = shape.Y;  % Latitude coordinates

[X,Y] = meshgrid(x, y);
[Latitude, Longitude] = ps2ll(X,Y); 

% Save as ISSM contour format for direct meshing
contour = [x(:), y(:)];

% Plot to verify
figure; plot(contour(:,1), contour(:,2), '-k');
axis equal;
title('Shapefile Contour');




filename = 'Ross_Ice_shelf.exp';
fid = fopen(filename, 'w');
fprintf(fid, 'iceshelf\n'); % Name of the domain
fprintf(fid, '%d\n', length(x)); % Number of points

for i = 1:length(x)
    fprintf(fid, '%f %f\n', x(i), y(i));
end

fclose(fid);

%%
% Initialize ISSM model
md = model();

% Mesh using BAMG (instead of Triangle)
hmin = 500;  % Minimum element size (500 m)
hmax = 5000; % Maximum element size (5 km)

md = bamg(model, 'domain', 'Ross_Ice_Shelf.exp', 'hmin', hmin, 'hmax', hmax);

% Plot mesh
plotmodel(md, 'data', md.mesh.x);
title('Generated Mesh using BAMG');


%%
