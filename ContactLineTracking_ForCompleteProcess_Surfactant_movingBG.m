


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code is designed for contact line tracking
% For surfactant, a switch criterion is first found and then used for track the inner film CL
% by Jing Shi, 29/10/2018

%%  PART ONE
% initiating parameters and variables

clear all;
%folder='H:\A-Project\A-Experiments\JingEXP\20170609-Et0.5-Si0.5Per-PS0.01Per@@\VideoFrames_II-4-side-withBg';
folder='D:\Jingdesktop\Durham Project\C3_C14E6(5mM)SDS(0.2mM)_4000fps_30v';
%%% go into the iamge folder
oldFolder=cd(folder);
% provide your video name (the part before numbers) and format
filetitle='C3_C14E6(5mM)SDS(0.2mM)_4000fps_30v';
Fileformat ='avi';
filename = [filetitle,'.',Fileformat];
frames = VideoReader(filename);
numFrames = get(frames, 'NumberOfFrames');

background_file_start=1; % with flat surface before drop landing 
background_file_end=5087;% the frame where inner film disappears

drying_start_frame=3; % whole drop landed
drying_end_frame=5057; % end of drying

inverse_frame_start=5053;% The frame where there is a very small inner ring
inverse_frame_end=5;% end frame of your interest


MovingBg=1; % whether to use a moving background; 0: no use
BG_Switchfile=1800; % the frame to switch between start background and moving/end background
                    % this frame should be at lease pass the spreading stage 
Bg_gap=25; % the gap between calculationg frame and the background frame ahead 



RingMask=1; % decide whether a ring mask is needed, 0 no need
frame_RingForm=5053;% frame used to decide the ring size
RingFrameRemove=4500; % the frame above which a mask is added to remove outsider rings
RingWidth=45; % the approximate ring width
skip_frames=0; % take frame+skip_frames+1, no skip: 0



%%% resolutions
resolution=0.40;%um/pixel
frame_rate=4000/(skip_frames+1);%fps 125 for test_01; water case 100fps, 50ethanol/water: \test_image_1: 3000fps. it has to be carefully done as camera could skip save frames
frame_dt=1/frame_rate;%time between neighbor frames /s

t_dry=(drying_end_frame-drying_start_frame)*frame_dt*1000; % total drying time in ms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% deposit shapes
%
background_sub_start=1; % case 0: no substraction, case 1: substraction
background_sub_end=1; % case 0: no substraction, case 1: substraction
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

%% Get the deposit ring diameter
% Remove background 
    switch background_sub_start
        case 0
            background_start=0;
        case 1
             background_start=read(frames,background_file_start);
             background_start=background_start(:,:,1);
    end
    
    
    switch background_sub_end
        case 0
            background_end=0;
        case 1
             background_end=read(frames,background_file_end);
             background_end=background_end(:,:,1);
    end
    % Get ROI from deposit image 
    if RingMask==1
    I0= read(frames,frame_RingForm);
    I0=I0(:,:,1);
    figure (1);
    cla
    imagesc(I0);
     colormap gray;
     axis equal;
     axis tight;
     
    I0=double(I0-background_start)+double(background_start-I0);% Background subtraction
%     I=abs(double(I-background));
    I0 = imrotate(I0,rotation_angle,'crop');
     I0=medfilt2(I0, [5 5]); % filter the noise
     I0= imgaussfilt(I0, 0.5); % Gaussian filter, change the second parameter if needed.
     I0=mat2gray(I0); % scale to 0-1
     figure(2)
     cla
     imagesc(I0);
     colormap gray;
     axis equal;
     axis tight;
     xlabel('x /pixel');
     ylabel('y /pixel');

%%% Edge detection of the deposit ring
           
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
     I =bwareaopen(I, 30);% remove smaller objects
     size_I=size(I);
     figure (3)
%    imshow(BW)
     imagesc(I);
     colormap gray;
     axis equal;
     axis([0 size_I(2) 0 size_I(1)]);
     xlabel('x /pixel');
     ylabel('y /pixel');
     
 %  determine the contour of the deposit ring,get the centroids of deposit ring and the diameter for analysis
     [n,m]=find(I); % in matrix row 'n' is coordinates' y and column 'm' is x           
     pklist=[m n];
     indCL=convhull(pklist(:,1),pklist(:,2)); 
     CL=pklist(indCL,:);

    % Plot the convex hull of each frame
    figure (1)
     hold on
     %plot(CL(:,1),CL(:,2),'ro','Markersize',2,'Linewidth',1);
     text(1,20,'Deposit Ring','color','r','Fontsize',10);
     
        %Fit a circle and extract the centre
        circ_fit=[CL ones(length(CL),1)]\[-(CL(:,1).^2+CL(:,2).^2)];
        xc_depositRing = -.5*circ_fit(1);
        yc_depositRing = -.5*circ_fit(2);
        R_c_depositRing  =  sqrt((circ_fit(1)^2+circ_fit(2)^2)/4-circ_fit(3));
        th=0:0.01:2*pi; %%%%Angle
        xy_c=[xc_depositRing+R_c_depositRing*cos(th) ; yc_depositRing+R_c_depositRing*sin(th)]';
        hold on
        plot(xy_c(:,1),xy_c(:,2),'-r','LineWidth',1)
        plot(xc_depositRing,yc_depositRing,'xr','Markersize',10,'Linewidth',1)
    end  
   % close all;
%% Get contact line of each frame (within R_c_depositRing)
frameNum=floor((inverse_frame_start-inverse_frame_end)/(skip_frames+1))+1;
frame_i=zeros(frameNum,1);
time=zeros(frameNum,1);
X_c=zeros (frameNum,1);
Y_c=zeros (frameNum,1);
R_C=zeros (frameNum,1);
    
count=0;
for i=inverse_frame_start:-(skip_frames+1):inverse_frame_end
    count=count+1;
    % read original image
    I0=read(frames,i);
    I0=I0(:,:,1);
    figure (1);
    cla
    imagesc(I0);
     colormap gray;
     axis equal;
     axis tight;
     xlabel('x /pixel');
     ylabel('y /pixel');
    if  MovingBg==1 && i>BG_Switchfile && i+Bg_gap < background_file_end 
    
        % read the image Bg_gap frame ahead as backgroud for substration
        I02=read(frames,i+Bg_gap);
        I02=I02(:,:,1);
        figure (2);
        cla
        imagesc(I02);
         colormap gray;
         axis equal;
         axis tight;
         xlabel('x /pixel');
         ylabel('y /pixel');

         % remove background, rotatoin
    %     if i>=frame_RingForm
    %         I0=double(I0-background_end)+double(background_end-I0);% Background subtraction
    %     elseif i> background_file_middle
    %         I0=double(I0-background_start)+double(background_start-I0);% Background subtraction
    %     end
            if i+Bg_gap < background_file_end
                I0=double(I0-I02)+double(I02-I0);        
            else 
                I0=double(I0-background_end)+double(background_end-I0);
            end
        
    else
         if i>=frame_RingForm
            I0=double(I0-background_end)+double(background_end-I0);% Background subtraction
         else
            I0=double(I0-background_start)+double(background_start-I0);% Background subtraction
         end
 
    end
     
    I0 = imrotate(I0,rotation_angle,'crop');
    I0=medfilt2(I0, [5 5]); % filter the noise
    I0= imgaussfilt(I0, 0.5); % Gaussian filter, change the second parameter if needed.
    I0=mat2gray(I0); % scale to 0-1
     figure(3)
     cla
     imagesc(I0);
     colormap gray;
     axis equal;
     axis tight;
     xlabel('x /pixel');
     ylabel('y /pixel');
     
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
 %    I = imbinarize(I0,'adaptive','ForegroundPolarity','bright');
     
     if background_sub_start==0
         I=imcomplement(I);
     end
     
     
     if i>=frame_RingForm
     I =bwareaopen(I,300);% remove smaller objects 
     else
     I =bwareaopen(I, 600);% remove smaller objects; larger threthold is used as the ring is large
     end
     
     I= imfill(I,'holes'); % fill holes
%      figure(3)
%      cla
%      imagesc(I);
%      colormap gray;
%      axis equal;
%      axis tight;
%      xlabel('x /pixel');
%      ylabel('y /pixel');
     
     % Create a circlur mask to remove the deposit ring
      if RingMask==1 && i>=RingFrameRemove
         [xx,yy] = ndgrid((1:size_I(1))-yc_depositRing,(1:size_I(2))-xc_depositRing);
         mask = double((xx.^2 + yy.^2)<(R_c_depositRing-RingWidth)^2); % creat a mask which is slightly smaller than the deposit ring
         I=I.*mask;
      end
     figure(4)
     cla
     imagesc(I);
     colormap gray;
     axis equal;
     axis tight;
     xlabel('x /pixel');
     ylabel('y /pixel');
     
 % Get contour of contact lines

if  MovingBg==1 && i>BG_Switchfile
    props = regionprops('table',I,'Centroid','MajorAxisLength','MinorAxisLength', 'ConvexHull');
    diameters = mean([props.MajorAxisLength props.MinorAxisLength],2);
    radii = diameters/2;
    [radii,Index]=sort(radii,'descend');
    props=props(Index,:);

  % Plot the convex hull of the biggest object
     CL=props.ConvexHull{1};
else
     [n,m]=find(I); % in matrix row 'n' is coordinates' y and column 'm' is x           
     pklist=[m n];
     indCL=convhull(pklist(:,1),pklist(:,2)); 
     CL=pklist(indCL,:);  
end

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
        Printfilename= ['_Frame',num2str(i)];
        print(fig,Printfilename,'-dtiff','-r100' );
        end 
        %pause ;

end
 time_dimentionless=1/t_dry.*time;
 Dc_real=2*R_C*resolution; % in um
 Dc_dimentionless=1/R_C(end-2).*R_C; % use R_C (end-2) to remove the influence of contact
 contactInfor=[frame_i,time,time_dimentionless,X_c,Y_c,R_C,Dc_real,Dc_dimentionless];
 % find local spikes of Diamter which is wrong data and delete them
%  [pks,indPeaks] = findpeaks(contactInfor(:,6));
%  contactInfor (indPeaks, :)=[];
%% Plot and print figures
% D-t
figure (5)
plot(time, Dc_real,'-o')
% set(gca, 'XLim', [0, 2000],'YLim', [0, 250])
xlabel('time / ms');
ylabel('Diameter /\mum');

fig = gcf;
fig.Color = 'white'; % set the background figure color
fig.InvertHardcopy = 'off';
iptsetpref('ImshowBorder','tight'); % Figures without any borders 
print(fig,'_D-t','-dtiff','-r100' );

% D-t (D-dimentionless)
figure (6)
loglog(time,Dc_dimentionless,'-o')
set(gca, 'XLim', [1, 2000],'YLim', [0.1, 2])
xlabel('time / ms');
ylabel('D/D_{0}');

fig = gcf;
fig.Color = 'white'; % set the background figure color
fig.InvertHardcopy = 'off';
iptsetpref('ImshowBorder','tight'); % Figures without any borders 
print(fig,'_D_D0-t','-dtiff','-r100' );


% D-t (D-dimentionless, t-dimentionless)
figure (8)
loglog(time_dimentionless,Dc_dimentionless,'-o')
set(gca, 'XLim', [0.001, 2],'YLim', [0.1, 2])
xlabel('t/t_{dry}');
ylabel('D/D_{0}');

fig = gcf;
fig.Color = 'white'; % set the background figure color
fig.InvertHardcopy = 'off';
iptsetpref('ImshowBorder','tight'); % Figures without any borders 
print(fig,'_D_D0-t_t_dry','-dtiff','-r100' );

% spreading part
figure (7)
loglog(time,Dc_dimentionless,'-o')
set(gca, 'XLim', [1, 100],'YLim', [1,1.8])
xlabel('time / ms');
ylabel('D/D_{0}');

fig = gcf;
fig.Color = 'white'; % set the background figure color
fig.InvertHardcopy = 'off';
iptsetpref('ImshowBorder','tight'); % Figures without any borders 
print(fig,'_D_D0-t_100ms','-dtiff','-r100' );

% D-t (D-dimentionless, t-dimentionless, semlog)
figure (9)
semilogy(time_dimentionless,Dc_dimentionless,'-o')
set(gca, 'XLim', [0.001, 1],'YLim', [0.1, 2])
xlabel('t/t_{dry}');
ylabel('D/D_{0}');

fig = gcf;
fig.Color = 'white'; % set the background figure color
fig.InvertHardcopy = 'off';
iptsetpref('ImshowBorder','tight'); % Figures without any borders 
print(fig,'_D_D0-t_t_dry','-dtiff','-r100' );

%% write to .csv file
save _contactInfor.mat contactInfor
str = date; % get the current date 
fileSavename = ['_contactInfor_',str,'.','csv'];
csvwrite(fileSavename,contactInfor);
%%
cd(oldFolder);
