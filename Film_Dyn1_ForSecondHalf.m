

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code is designed for fringe analysis 
% by Lisong Yang, 02-02-2016 for fringe dynamic analysis

%%  PART ONE
% initiating parameters and variables

clear all;
folder='H:\M-Matlab codes\Film interferometry\Surfactants';
drop_fringes_file_1025 % load file information
%%
background_sub=1; % case 0: no substraction, case 1: substraction
fringe_enhance=0; % background smooth. case 0: no, case 1:; case 2:
intensity_threshold_contact=50; %*!Line detection. Please use cross section test to set this value

bright_ring_threshold=1;% identifying bright rings
dark_ring_threshold=1; % identifying dark rings

smooth_points_t=3; % for tracing over time
threshold_peaks_t=10;  %set the criteria to distinguish peaks and valleys

smooth_points_r=3; %
threshold_peaks_r=10;  %set the criteria to distinguish peaks and valleys

ROI_auto=0; % case 0: auto determination, case 1: manual determination
contact_fixed=0;
rotation_angle=0;% %2.7 for square; %deg

oldFolder=cd(folder);

%% Remove background 
    switch background_sub
        case 0
            background=0;
        case 1
             filename = [filetitle,num2str(background_file_start),'.',fileformat];
             filename=fullfile (folder, filename);
             background = imread(filename); 
             background=background(:,:,1);
    end
    % Get ROI from deposit image
    filename = [filetitle,num2str(inverse_frame_start),'.',fileformat];
    I0= imread(filename);
    I0=I0(:,:,1);
    figure (1);
    cla
    imagesc(I0);
     colormap gray;
     axis equal;
     axis tight;
     
    I0=double(I0-background)+double(background-I0);% Background subtraction
%     I=abs(double(I-background));
    I0 = imrotate(I0,rotation_angle,'crop');
     I0=medfilt2(I0, [5 5]); % filter the noise
     I0=mat2gray(I0); % scale to 0-1
     figure(2)
     cla
     imagesc(I0);
     colormap gray;
     axis equal;
     axis tight;
     xlabel('x /pixel');
     ylabel('y /pixel');

%% Edge detection of the deposit ring
           
 % Edge dection  (the brigtest part)
%     Gx = [1 -1; 1 -1]; 
%     Gy = [1 1; -1 -1];
% 
%     edgepic1 = abs(conv2(I0,Gx,'same'));
%     edgepic2 = abs(conv2(I0,Gy,'same'));
% 
%     I= sqrt(edgepic1(1:end-1,1:end-1).^2+edgepic2(1:end-1,1:end-1).^2);  
    
% use Matlab build-in edge detection algrithm
    %I=edge(I0,'canny');                 
     I=edge(I0,'sobel'); 
     I =bwareaopen(I, 10);% remove smaller objects
     size_I=size(I);
     figure (3)
%    imshow(BW)
     imagesc(I);
     colormap gray;
     axis equal;
     axis([0 size_I(2) 0 size_I(1)]);
     xlabel('x /pixel');
     ylabel('y /pixel');
     
 %%   determine the contour of the deposit ring,get the centroids of deposit ring and the diameter for analysis
     [n,m]=find(I); % in matrix row 'n' is coordinates' y and column 'm' is x           
     pklist=[m n];
     indCL=convhull(pklist(:,1),pklist(:,2)); 
     CL=pklist(indCL,:);

    % Plot the convex hull of each frame
     hold on
     plot(CL(:,1),CL(:,2),'ro','Markersize',2,'Linewidth',1);
     text(1,20,'Deposit Ring','color','white','Fontsize',10);
     
        %Fit a circle and extract the centre
        circ_fit=[CL ones(length(CL),1)]\[-(CL(:,1).^2+CL(:,2).^2)];
        xc_depositRing = -.5*circ_fit(1);
        yc_depositRing = -.5*circ_fit(2);
        R_c_depositRing  =  sqrt((circ_fit(1)^2+circ_fit(2)^2)/4-circ_fit(3));
        th=0:0.01:2*pi; %%%%Angle
        xy_c=[xc_depositRing+R_c_depositRing*cos(th) ; yc_depositRing+R_c_depositRing*sin(th)]';
        hold on
        plot(xy_c(:,1),xy_c(:,2),'-g','LineWidth',2)
        plot(xc_depositRing,yc_depositRing,'xg','Markersize',10,'Linewidth',2)
      
    close all;
%% Get contact line of each frame (within R_c_depositRing)
frameNum=floor((inverse_frame_start-inverse_frame_end)/(skip_frames+1))+1;
frame_i=zeros(frameNum,1);
time=zeros(frameNum,1);
X_c=zeros (frameNum,1);
Y_c=zeros (frameNum,1);
R_C=zeros (frameNum,1);
background_sub_end=1;
    switch background_sub_end
        case 0
            background_end=0;
        case 1
             filename = [filetitle,num2str(background_file_end),'.',fileformat];
             filename=fullfile (folder, filename);
             background_end = imread(filename); 
             background_end=background_end(:,:,1);
    end
    
count=0;
for i=inverse_frame_start:-(skip_frames+1):inverse_frame_end
    count=count+1;
    % read original image
    filename = [filetitle,num2str(i),'.',fileformat];
    I0= imread(filename);
    I0=I0(:,:,1);
    figure (1);
    cla
    imagesc(I0);
     colormap gray;
     axis equal;
     axis tight;
     xlabel('x /pixel');
     ylabel('y /pixel');
    
    % remove background, rotatoin
    I0=double(I0-background_end)+double(background_end-I0);% Background subtraction
    I0 = imrotate(I0,rotation_angle,'crop');
    I0=medfilt2(I0, [5 5]); % filter the noise
    I0=mat2gray(I0); % scale to 0-1
%      figure(2)
%      cla
%      imagesc(I0);
%      colormap gray;
%      axis equal;
%      axis tight;
%      xlabel('x /pixel');
%      ylabel('y /pixel');
     
% % Method1:Edge dection  (the brigtest part)
%     Gx = [1 -1; 1 -1]; 
%     Gy = [1 1; -1 -1];
% 
%     edgepic1 = abs(conv2(I0,Gx,'same'));
%     edgepic2 = abs(conv2(I0,Gy,'same'));
% 
%     I= sqrt(edgepic1(1:end-1,1:end-1).^2+edgepic2(1:end-1,1:end-1).^2);  
    
% % Method1: use Matlab build-in edge detection algrithm
%     %I=edge(I0,'canny');                 
%      I=edge(I0,'sobel'); 
%      I =bwareaopen(I, 20);% remove smaller objects
%      figure (3)
%      imagesc(I);
%      colormap gray;
%      axis equal;
%      axis([0 size_I(2) 0 size_I(1)]);
%      xlabel('x /pixel');
%      ylabel('y /pixel');

% Method 3:  globle threshold method
     level = graythresh(I0);
     I = im2bw(I0, level);
     if background_sub==0
         I=imcomplement(I);
     end
     I =bwareaopen(I, 50);% remove smaller objects
%      figure(3)
%      cla
%      imagesc(I);
%      colormap gray;
%      axis equal;
%      axis tight;
%      xlabel('x /pixel');
%      ylabel('y /pixel');
     
     % Create a circlur mask to remove the deposit ring
     [xx,yy] = ndgrid((1:size_I(1))-yc_depositRing,(1:size_I(2))-xc_depositRing);
     mask = double((xx.^2 + yy.^2)<(R_c_depositRing-5)^2); % creat a mask which is slightly smaller than the deposit ring
     I=I.*mask;
     figure(4)
     cla
     imagesc(I);
     colormap gray;
     axis equal;
     axis tight;
     xlabel('x /pixel');
     ylabel('y /pixel');
     
     % Get contour of contact lines
     [n,m]=find(I); % in matrix row 'n' is coordinates' y and column 'm' is x           
     pklist=[m n];
     indCL=convhull(pklist(:,1),pklist(:,2)); 
     CL=pklist(indCL,:);

    % Plot the convex hull of each frame
     figure (1)
     hold on
     %plot(CL(:,1),CL(:,2),'ro','Markersize',2,'Linewidth',1);
     message=sprintf('Frame%4.0f',i);
     text(10,20,message,'color','r');
          
        %Fit a circle and extract the centre
        circ_fit=[CL ones(length(CL),1)]\[-(CL(:,1).^2+CL(:,2).^2)];
        xc_c = -.5*circ_fit(1);
        yc_c = -.5*circ_fit(2);
        R_c_c  =  sqrt((circ_fit(1)^2+circ_fit(2)^2)/4-circ_fit(3));
        th=0:0.01:2*pi; %%%%Angle
        xy_c=[xc_c+R_c_c*cos(th) ; yc_c+R_c_c*sin(th)]';
        
        plot(xy_c(:,1),xy_c(:,2),'-r')
        plot(xc_c,yc_c,'xr')
        frame_i(count)=i;
        time(count)=frame_dt*(i-drying_start_frame)*1000;% in ms
        X_c(count)=xc_c;
        Y_c(count)=yc_c;
        R_C(count)= R_c_c;
          
          %%% display the processed percentage and print iamges
        if mod(i-inverse_frame_start,100)==0
        message=sprintf('processed =%5.1f %%',abs((i-inverse_frame_start))/frameNum*100);
        disp(message);
        % print images
        fig = gcf;
        fig.Color = 'white'; % set the background figure color
        fig.InvertHardcopy = 'off';
        iptsetpref('ImshowBorder','tight'); % Figures without any borders 
        Printfilename= ['_',filetitle,num2str(i)];
        print(fig,Printfilename,'-dtiff','-r100' );
        end 
        %pause (0.1);

end
 time_dimentionless=1/t_dry.*time;
 Dc_real=2*R_C*resolution; % in um
 contactInfor=[frame_i,time,time_dimentionless,X_c,Y_c,R_C,Dc_real];
%% Plot and print figures
% D-t
figure (5)
plot(time, Dc_real,'-o')
xlabel('time / ms');
ylabel('Diameter /pL');

fig = gcf;
fig.Color = 'white'; % set the background figure color
fig.InvertHardcopy = 'off';
iptsetpref('ImshowBorder','tight'); % Figures without any borders 

print(fig,'_D-time','-dtiff','-r100' );

 
%% write to .csv file
str = date; % get the current date 
fileSavename = ['_contactInfor_',str,'.','csv'];
csvwrite(fileSavename,contactInfor);
%%
cd(oldFolder);
