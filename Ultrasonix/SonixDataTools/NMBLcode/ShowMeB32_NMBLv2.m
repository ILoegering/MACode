% Clearing the environment
clc
clear
if ishandle(1010), close 1010, end
format compact

% Making the other functions accessible to the environment
% addpath('/Volumes/External/NMBL/TendonSurfaceArea/Ultrasonix/SonixDataTools/DataReaders/RPread/');
% addpath('/Volumes/External/NMBL/TendonSurfaceArea/Ultrasonix/SonixDataTools/ImageProcessing/');
% addpath('/Volumes/External/NMBL/TendonSurfaceArea/Ultrasonix/SonixDataTools/PlotFunctions/');
% addpath('/Volumes/External/NMBL/TendonSurfaceArea/Ultrasonix/SonixDataTools/Misc/');

%% Reading the data
dpath = 'D:\UW_NMBL\vibrationInVivo\AchillesSpatialVariation\Achilles_CSA\';      % File path
subject = 'achilles';                               % Subject
trial = '0mm';                                     % Image position
% folder = strcat(subject,'\');                       % Subject's folder
folder = '';
filename = strcat(subject,'_',trial);
filetype = '.b32';                                  % File type
[Data, Header] = RPread(strcat(dpath,folder,filename,filetype), '6.0.3');

%% Plotting the ultrasound .b32 image

% Creating a new plot
hF = figure(1010);
hA = axes;

% Plotting the first frame
% plot_SonixRP(Data(:, :, 1), Header , [hA -1], 1);
plot_SonixRP(Data(:, :, 1), Header , [hA hA], 1);
hold on

%% Selecting and plotting the points for the outline of the tendon
title('Outline cross-sectional area')
disp('Outline cross-sectional area')
n = 40;                 % number of points used in outline
for i = 1 : n
    try
        [nx(i,:),ny(i,:)] = ginput(1);    % Coordinates of outline in pixels
        plot(nx(i),ny(i),'rx')
    catch                                % Use Ctrl+c to break
        nx(i) = nx(1); ny(i) = ny(1);     % Closing polygon
        break
    end
end

% Plotting the outline of the points selected
plot(nx,ny,'r','linewidth',1.25);
plot(nx,ny,'or','markersize',3);
conv = 87/1000;         % Image is 87 microns per pixel (find this in html file also saved with .b32)
x = conv.*nx;
y = conv.*ny;
coords = [x y];         % Coordinates converted from pixels to mm

% Calculating area
[XX,YY] = meshgrid(1:size(Data,1),1:size(Data,2));
[inOutline,onOutline] = inpolygon(XX,YY,nx,ny);
conv = 87/1000;         % Image is 62 microns per pixel
tendonArea = sum(sum(inOutline))*conv^2;

%% Selecting and plotting the points for the height and width of tendon
title('Choose two points for height')
disp('Choose two points for height')
for i = 1 : 2
    [hx(i,:),hy(i,:)] = ginput(1);    % Coodrinates of the tendon height
    plot(hx(i),hy(i),'yx')
end
title('Choose two points for width')
disp('Choose two points for width')
for i = 1 : 2
    [wx(i,:),wy(i,:)] = ginput(1);    % Coordinates of the tendon width
    plot(wx(i),wy(i),'yx')
end

% height
plot(hx,hy,'--y','linewidth',0.75);
plot(hx,hy,'*y','markersize',5);
heightx = conv.*hx;
heighty = conv.*hy;
h = [heightx heighty];  % Coordinates converted from pixels to mm
height = pdist(h,'euclidean');

% width
plot(wx,wy,'--y','linewidth',0.75);
plot(wx,wy,'*y','markersize',5);
widthx = conv.*wx;
widthy = conv.*wy;
w = [widthx widthy];    % Coordinates converted from pixels to mm
width = pdist(w,'euclidean');

%% Calculating the 2D area of tendon outline 
%     % For each line segment, this formula makes a rectangle and sums the 
%     % area in a counterclockwise manner. The width is the distance from the 
%     % midpoint of the x verticies to the y axis of the image and the height 
%     % is the difference of the y verticies. When another segment comes
%     % across the rectangle, the excess area of the first rectangle is
%     % subtracted and the area is only of the inside of the polygon.
%     % http://www.mathopenref.com/coordpolygonarea2.html
% 
% area = 0;               % Accumulates area in the loop
% numPoints = length(x);  % The number of verticies
% j = numPoints-1;        % The last vertex is the previous one to the first
% 
%   for i=1:numPoints
%     area = area +  (x(j)+x(i)) * (y(j)-y(i)); 
%       j = i;
%   end
%   
% tendonArea = abs(area/2);   % Divided by 2 to use the midpoint of the distances

%% Saving image with area values
title('Select locations for labels')
disp('Select locations for labels')
gtext(strcat('Achilles Area = ','  ',num2str(tendonArea),' mm^2'),...
    'Color','white','FontSize',8);
gtext(strcat('Height = ',' ',num2str(height),' mm'),'Color','white',...
    'FontSize',8);
gtext(strcat('Width = ',' ',num2str(width),' mm'),'Color','white',...
    'FontSize',8);
% cd /Volumes/External/NMBL/TendonSurfaceArea/Images/Processed/  % Changing directory 
cd D:\UW_NMBL\vibrationInVivo\AchillesSpatialVariation\Achilles_CSA\processed\  % Changing directory 

% Formatting image to print
H=7;W=11;
set(hF,'PaperUnits','inches')
set(hF,'PaperOrientation','portrait')
set(hF,'PaperSize',[H,W])
set(hF,'PaperPosition',[0,0,W,H])
% print(hF,'-dpng','AASY0001_acc1.jpg','-r600','-S3600,s2600','-tight');
print(hF,'-dpng',strcat(filename,'.jpg'));
% cd /Volumes/External/NMBL/TendonSurfaceArea/Ultrasonix/SonixDataTools/ExampleScripts/
cd D:\UW_NMBL\vibrationInVivo\

% close 1010