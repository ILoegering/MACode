function[tendonCoords] = tendonDepth(Data,Header,xConv,yConv,tw,framesOfInterest)
%tendonDepth - get coordinates of tendon mid-line in ultrasound images
%Superficial and deep edges of tendon are selected manually, and the
%coordinates of the mid-line of the tendon are calculated.
%
%
% Syntax:  [tendonCoords] = tendonDepth(Data,Header,xConv,yConv,tw,framesOfInterest)
%
% Inputs:
%    Data (required) - double matrix (numWPix x numDPix x numFrames)
%           ultrasound data; first output of Ultrasonix function RPread.m
%    Header (required) - struct (1 x 1)
%           ultrasound header info; second output of Ultrasonix function
%           RPread.m
%    xConv (required) - double (1 x 1)
%           conversion factor from pixels to length in meters/pixel for
%           depth dimension (x-axis)
%    yConv (required) - double
%           conversion factor from pixels to length in meters/pixel for
%           width dimension (y-axis)
%    tw (required) - double (1 x 1)
%           transducer width (along the tendon) in meters
%    framesOfInterest (required) - double array (1 x numFramesOfInterest)
%           frames of ultrasound data for which tendon coordinates are to
%           be determined
%                       
%
% Outputs:
%    tendonCoords - struct (numFramesOfInterest x 1)
%           tendon coordinates in ultrasound reference frame
%
% Other m-files required: plot_SonixRP.m
% Subfunctions: none
% MAT-files required: none
%
% Author: Jack Martin
% Editor: Isaac Loegering
% UW Neuromuscular Biomechanics Lab
% University of Wisconsin-Madison
% 1513 University Ave, Rm 3046
% Madison, WI 53706
% email: jack.martin1313@gmail.com
%        isaacloegering@gmail.com
% A long time ago in a biomechanics lab far, far away; Last revision: 21-May-2019
%------------- BEGIN CODE --------------
figure
set(gcf,'WindowState','Maximize');
hA = axes;
f = 0;
% For each frame of interest
for frame = framesOfInterest
    f = f + 1;
    
    % Get image dimensions
    dim = size(Data(:, :, frame));
    
    % Loop until satisfied with resulting tendon mid-line
    while (1)
        % Free memory
        clearvars -except frame f hA hF Data Header tendonCoords xConv yConv tw dim
        
        % Plot ultrasound .b32 image
        plot_SonixRP(Data(:, :, frame), Header , [hA hA], 1);
        
        % Set y-axis to display normally (+y-direction up)
        set(gca,'YDir','normal')
        hold on
        
        % Selecting and plotting the points for the tendon boundaries
        % Superficial Boundary
        title({'Select points spanning superficial tendon boundary', 'click outside image when finished'})
        fprintf('\n\nSelect 7 points spanning superficial tendon boundary \n click outside image when finished')
        i = 1;
        while (1)
            % Get point
            [x,y] = ginput(1);
            % If point is outside image, break
            if (x > dim(2) || x < 0 || y > dim(1) || y < 0)
                break
            end
            % If point is inside image, store it and loop
            nxU(i,:) = x; nyU(i,:) = y;    % Coordinates of outline in pixels
            plot(nxU(i),nyU(i),'rx','linewidth',2,'markersize',6)
            i = i + 1;
        end
        % Deep Boundary
        title({'Select points spanning deep tendon boundary', 'click outside image when finished'})
        fprintf('\n\nSelect 7 points spanning deep tendon boundary \n click outside image when finished')
        i = 1;
        while (1)
            % Get point
            [x,y] = ginput(1);
            % If point is outside image, break
            if (x > dim(2) || x < 0 || y > dim(1) || y < 0)
                break
            end
            % If point is inside image, store it and loop
            nxL(i,:) = x; nyL(i,:) = y;    % Coordinates of outline in pixels
            plot(nxL(i),nyL(i),'rx','linewidth',2,'markersize',6)
            i = i + 1;
        end
        
        % Plotting the outline of the points selected
        plot(nxU,nyU,'ro:',nxL,nyL,'ro:','markersize',6);
        xU = xConv.*nxU;
        yU = yConv.*nyU;
        xL = xConv.*nxL;
        yL = yConv.*nyL;
        
        % Interpolating tendon boundaries
        Y = linspace(0,tw,20);
        XU = spline(yU,xU,Y);
        XL = spline(yL,xL,Y);
        
        % Calculate tendon mid-line as average of superficial and deep
        % boundaries
        XC = mean([XU; XL],1);
        
        % Plot superficial and deep boundaries and tendon mid-line
        plot(XU/xConv,Y/yConv,'c-',XL/xConv,Y/yConv,'c-')
        plot(XC/xConv,Y/yConv,'c-','linewidth',3)
        hold off
        
        % Approve tendon-midline or repeat procedure
        title({'Left click if satisfied', 'right click to try again'})
        fprintf('\n\nLeft click if satisfied \n right click to try again\n')
        [~,~,clickVal] = ginput(1);
        if clickVal == 1
            break
        end
    end
    
    % Store coordinates of tendon mid-line
    tendonCoords.depth{f} = XC;
    tendonCoords.axialLocation{f} = Y;
    tendonCoords.frame(f) = frame;
end
close
%------------- END OF CODE --------------
end
