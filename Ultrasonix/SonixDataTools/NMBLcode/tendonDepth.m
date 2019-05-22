function[tendonCoords] = tendonDepth(directory,filename,framesOfInterest)

if ishandle(1010), close 1010, end

%% Reading the data

if nargin == 0
    [filename directory]=uigetfile(['*.b32']);
end
% if nargin == 0 || nargin == 2
%     decimation_factor = 10; % default decimation factor for downsampling
% end


[Data, Header] = RPread(strcat(directory,filename), 1000);

% Trim 50 px border
border_size = 50;
Data = Data(border_size+1:end-border_size,border_size+1:end-border_size,:);
% Flip y-axis direction (reflect data across x-axis so origin is in lower left corner when YDir is normal)
Data = flipud(Data);

%% Plotting the ultrasound .b32 image

hF = figure(1010);
hA = axes;

f = 0;

% for frame = 1 : decimation_factor : size(Data,3)
for frame = framesOfInterest
    f = f + 1;
    
    for attempt = 1 : 100
        clearvars -except frame f hA hF Data Header directory filename tendonCoords
        
        plot_SonixRP(Data(:, :, frame), Header , [hA hA], 1);
        % Set y-axis to display normally (+y-direction up)
        set(gca,'YDir','normal')
        hold on
        
        % Selecting and plotting the points for the tendon boundaries
        title({'Select points spanning superficial tendon boundary', 'press ctrl+c when finished'})
        fprintf('\n\nSelect 7 points spanning superficial tendon boundary \n press ctrl+c when finished')
        n = 40;                 % number of points used in outline
        for i = 1 : n
            try
                [nxU(i,:),nyU(i,:)] = ginput(1);    % Coordinates of outline in pixels
                plot(nxU(i),nyU(i),'rx','linewidth',2,'markersize',6)
            catch                                % Use Ctrl+c to break
                break
            end
        end
        
        title({'Select points spanning deep tendon boundary', 'press ctrl+c when finished'})
        fprintf('\n\nSelect 7 points spanning deep tendon boundary \n press ctrl+c when finished')
        for i = 1 : n
            try
                [nxL(i,:),nyL(i,:)] = ginput(1);    % Coordinates of outline in pixels
                plot(nxL(i),nyL(i),'rx','linewidth',2,'markersize',6)
            catch                                % Use Ctrl+c to break
                break
            end
        end
        
        % Plotting the outline of the points selected
        plot(nxU,nyU,'ro:',nxL,nyL,'ro:','markersize',6);
        conv = 87/1000;         % Image is 87 microns per pixel (find this in html file also saved with .b32)
        xU = conv.*nxU;
        yU = conv.*nyU;
        xL = conv.*nxL;
        yL = conv.*nyL;
        
        % Interpolating tendon boundaries
        Y = linspace(0,46.98,20);
        XU = spline(yU,xU,Y);
        XL = spline(yL,xL,Y);
        
        XC = mean([XU; XL],1);
        
        plot(XU/conv,Y/conv,'c-',XL/conv,Y/conv,'c-')
        plot(XC/conv,Y/conv,'c-','linewidth',3)
        
        title({'Left click if satisfied', 'right click to try again'})
        fprintf('\n\nLeft click if satisfied \n right click to try again')
        
        [~,~,clickVal] = ginput(1);
        
        if clickVal == 1
            break
        end
        
        hold off
    end
    
    tendonCoords.depth{f} = XC
    tendonCoords.axialLocation{f} = Y
    tendonCoords.frame(f) = frame
end

