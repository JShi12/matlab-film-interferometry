
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code is designed for fringe analysis for droplet with non-fixed CLs
% The centre and CL are inputs; results from ContactLineTracking code
% Adapted by Jing Shi from fringes analysis code by Lisong, Teresa and Colin


%%%%%%%%%%%initiating parameters and variables
% RUN drop_fringes_file.m first - which is save in the folder for individual experiment
% 
%%%%%%%%% Some comments on figures and variables
% figure 1: original image (or edge detected in case of bank_edge=1), 
% ..........ginput, cross section close to the centre of the droplet 
% figure 2: the result of cross section of the image, I(x) to be used to
% .........define the variable 'intensity_threshold_contact' 
% figure 3: cropped image (in case of circle fringes, derive the centre, in green cross, and
% .......... radius of the droplet), the centre of the image is also plotted in red circle
% figure 4: showing the middle term how the centre of a circle is automatically detected 
% figure 5: the process of creating I_3D(x,y,t), x and y are in pixels, 
% ..........t is in frame. Note that the first frame in I_3D is the end of
% .......... drying. The centre of image (or in case of circle drop the centre of droplet and contact line) 
% .........is plotted on the last showed image 
% figure 6: L-t at the centre of drop/image (intensity) through the analysed frames (index)
% figure 7: h_c-t at the centre of drop (H=lamda/(2n))
% figure 8: L(t)-h_c(t) at the centre of drop. This shows how intensity of
% interference signal associates with the film thickness. The experimental
% data fitted with damped sinsoidal funtion x1+x2*exp(x3*h^2)sin(x4*h+x5)
% figure 9: the image that is near the start of drop drying with the centre
% showing intensity peak or valley only!!

% I - 2D image data, x*y array
% L - 1D data, with index of frame (step 1) - used for centre drop intensity changes
% I_3D -3D data, x*y*t array image data with different frame. frame starts from dried drop
% 

% maxtab and mintab - 2D data x (either position or time) at max/min value of y (intensity with fringes - bright or dark)
% h_c - converted film thickness at the centre of drop from maxtab and mintab. [frames, thickness]
% cross_r - [x y r intensity], cross section pixels
% profile - [r h] - peaks and valleys only with derived thickness
% H_t0 - initial profile at defined t0 [x y r I h]
% H_r_t - 3D array [t h] [r]

% fringes_t [frames at peaks/valleys; intensity]
% height_1 [frames intensity height1], traced from flat substrate
% height_2 [frames height2], traced from drop near the start of drying - coaxial rings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Provide files Information
clear all

folder='R:\Jing Shi\A2-Project\EXP\20181018_Spreading_C14E6\B4_0.5mM_4000fps_30v';
filetitle='B4_0.5mM_4000fps_30v'; %50fps
fileformat ='avi';

oldFolder=cd(folder);
filename = [filetitle,'.',fileformat];

frames = VideoReader(filename);%this part is new
numFrames = get(frames, 'NumberOfFrames');

% Load contact infor
load _contactInfor.mat;
FrameAnalysis=contactInfor(:,1); % Frame number 
xc_c=contactInfor(:,4); % x-coordinate of the centre
yc_c=contactInfor(:,5); % y -coordinate of the centre
R_c=contactInfor(:,6); % Radius of the CL
background_file=4772; % no feature or with only bank structure (OLED case) 

drying_start_frame=3; % make sure it is the first frame landing, one frame after the background
drying_end_frame=4772; % end of drying,  this to calculate the drying time

analysis_start=FrameAnalysis(1);% start frame of your interst; 
analysis_end=3200;%end frame of your interest
intialPeakBrt=1; % the bright peak counts before the analysis_start frame (change of bright and dark in the centre)
intialPeakDark=1; % the dark peak counts before the analysis_start frame (change of bright and dark in the centre)
skip_frames=3; % take frame+skip_frames+1, no skip: 0

%%%%%%%%%%%%%%%%%%%
%%% resolutions
frame_rate=4000;
resolution=0.40;%um/pixel
frame_rate=frame_rate/(skip_frames+1);%fps 125 for test_01; water case 100fps, 50ethanol/water: \test_image_1: 3000fps. it has to be carefully done as camera could skip save frames
frame_dt=1/frame_rate;%time between neighbor frames /s
drying_time_in_frame=drying_end_frame-drying_start_frame;
drying_time=(drying_end_frame-drying_start_frame)/frame_rate;%s


% whether initial imgage has clear edge/structure

rotation_angle=0;% deg

%optics
lambda=0.470; %um blue LED thorlab 470nm
n_index=1.41; % 1.41 for C14E6 (Zhmud and Tiberg, 2005)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% intensity threshold
background_sub=0;
gaussian=1; % gaussion filter

fringe_enhance=0; % background smooth. case 0: no, case 1:; case 2:


bright_ring_threshold=1;% identifying bright rings
dark_ring_threshold=1; % identifying dark rings

%%%% other impportant variables that are needed to adjust below (search the
%%%% lines that are related...
% y1=smooth(double(intensity),smooth_points,'lowess'); % the second parameter is important. Adjust it according to your image contrast and noise level
% [maxtab, mintab] = peakdet(intensity, threshold_peaks, d) % keep the second parameter as low as possible. but be aware of noises
smooth_points_t=6; % for tracing over time
threshold_peaks_t=6;  %set the criteria to distinguish peaks and valleys

smooth_points_r=6; % for tracing cross section
smooth_points_r_head=2;
smooth_points_r_tail=2;
threshold_peaks_r=4;  %set the criteria to distinguish peaks and valleys

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get the size of one of the frames. Note no crop here, the videos are supposed to be croped already
  background = read(frames,background_file);
  background=background(:,:,1);
  background = imrotate(background(:,:,1),rotation_angle,'crop');
  size_I=size(background);
%% run through frames from dried towards the start of landing of droplet 
% save into I_3D
frame_i=1;
frame_i_counterPart=1;
frames_count=floor((analysis_start-analysis_end)/(skip_frames+1))+1;
I_3D=zeros(size_I(1),size_I(2),frames_count);
xc_c2=zeros(frames_count,1);
yc_c2=zeros(frames_count,1);
R_c2=zeros(frames_count,1);
for drop_frame=analysis_start:-(skip_frames+1):analysis_end
 
    I = read(frames,drop_frame);
    I=I(:,:,1);
  
       if background_sub   
             I=double(background-I)+double(I-background);%%%%%%%%%%%%%%%%%% Background subtract
             %I=abs(double(I-background));%%%%%%%%%%%%%%%%%% Background subtract
       end         
 
       if gaussian
%           filt_gauss=[2 4 5 4 2 ; 4 9 12 9 4; 5 12 15 12 5; 4 9 12 9 4; 2 4 5 4 2]./159;
%           I=conv2(filt_gauss,I);
           I=imgaussfilt(I,0.5); % the second parameter is Standard deviation of the Gaussian distribution. 0.5 is the default value in matlab
       end
          
     % rotation and crop    
    I = imrotate(I,rotation_angle,'crop'); %%% rotation
%     if gaussian
%       I=I(lefty_min+3:righty_max+3,leftx_min+3:rightx_max+3); %%% cropping
%     else
%       I=I(lefty_min:righty_max,leftx_min:rightx_max); %%% cropping
%     end
%     size_I=size(I)  
   
    I_3D(:,:,frame_i)=I;  % This contains information of frames in reverse order, i.e., I (:,:,1) for the drying end frame
    
    figure(5)
    cla
    imagesc(I_3D(:,:,frame_i));
    colormap gray; 
    axis equal;
    axis([0 size_I(2) 0 size_I(1)]);
  
     message=sprintf('Frames:%4.0f; Frame reversed:%4.0f; Time: %4.2f ms',drop_frame,frame_i,-(frame_i-1)*frame_dt*1000);
     text(10,20,message,'color','yellow');

    if background_sub
        text(10,40,'Orginal image background off','color','white');
    else
        text(10,40,'Orginal image','color','white');
    end
   % Plot the CL and the centre; 
   % Comment this section to save time
    hold on
  
    th=0:0.01:2*pi; %%%%Angle
    xy_c=[xc_c(frame_i_counterPart)+R_c(frame_i_counterPart)*cos(th) ; yc_c(frame_i_counterPart)+R_c(frame_i_counterPart)*sin(th)]';

    plot(xy_c(:,1),xy_c(:,2),'-r');
    plot(xc_c(frame_i_counterPart),yc_c(frame_i_counterPart),'xr');
    
   %%% Save the centres and radii of used frames
    xc_c2(frame_i)=xc_c(frame_i_counterPart);
    yc_c2(frame_i)=yc_c(frame_i_counterPart);
    R_c2(frame_i)=R_c(frame_i_counterPart);
    frame_i=frame_i+1;
    frame_i_counterPart=frame_i_counterPart+skip_frames+1;
    
end

%% Get the centre height by tracing over time (start from the end of drying to the early stage)
% the centre of the droplet is defined by the cropped image or in the case
% of circle droplet by the fitting, xc_c, yc_c
% converted film thickness vs time at the drop centre is saved into h_c
%!! This is for cases where the height change in one diection
cd(oldFolder);
L=zeros(frames_count,1); 
maxtab=[];
mintab=[];
h_c=[];
frame_t=1:frames_count;

%  test an arbitary point
%          figure(5)
%          [px,py,button]=ginput(1);
%  
%          for k=1:frames_count             
%                L=[L;I_3D(round(py),round(px),k)];     
%          end

        for k=1:frames_count             
               L(k)=I_3D(round(yc_c2(k)),round(xc_c2(k)),k); % The intensity value of the centre
        end
%         
        figure(6)
        cla;
         
        L1=smooth(L, smooth_points_t,'lowess'); % second variable is the smooth point, it is test here and used in below
       
               plot(L,'-b');
               hold on
               plot(L1,'r');
               xlabel('t /Frame');
               ylabel('Intensity');
               set(gca,'XMinorTick','on');
               axis tight;
                        
               [maxtab, mintab] = peakdet(L1, threshold_peaks_t, frame_t);
               % maxtab obtains (frame_t, peak intensity) and mintab obtains (frame_t, valley intensity)
             
               plot(mintab(:,1), mintab(:,2), 'g*');
               plot(maxtab(:,1), maxtab(:,2), 'r*'); 
             
           if ~isempty(mintab) & ~isempty(maxtab)
             for i=1:length(maxtab)
               h_c=[h_c; maxtab(i,1) (i+intialPeakBrt-1)*lambda/2/n_index];             
             end
             
             for i=1:length(mintab)
               h_c=[h_c; mintab(i,1) (i+intialPeakDark-1/2)*lambda/2/n_index];             
             end
             
             h_c=sort(h_c); % 2D data [frames, height/um that is related with peaks and valleys]
             
             I_c=[h_c(:,2) L(round(h_c(:,1)))]; % 2D data [thickness intensity]
             
             figure(7)
             cla
             plot( h_c(:,1),  h_c(:,2), 'g*');
             xlabel('t /Frame');
             ylabel('Thickness /\mum');
             axis tight;
             
%              % lsqcurvefit I-thickness
%              xdata=I_c(:,1);
%              ydata=I_c(:,2);
%              fun = @(x,xdata) x(1)+x(2)*exp(x(3)*xdata.^2).*sin(x(4)*xdata+x(5));
%              x0=[mean(ydata),(max(ydata)-min(ydata))/2,-0.8,2*pi/0.15,pi/2]; % initial parameters, estimated from original data
%              % options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt');
%              lb = [mean(ydata)-2,80,-0.2,40,1];
%              ub = [mean(ydata)+2,88,-0.14,45,2];
%              x = lsqcurvefit(fun,x0,xdata,ydata,lb,ub);
%              thickness=0:0.01:4;
% 
%              figure(8)
%              cla
%              plot(I_c(:,1),I_c(:,2), 'g*');
%                xlabel('Thickness /\mum');
%                ylabel('Intensity');
%                hold on
%                plot(thickness,fun(x,thickness),'r-');    
%                axis tight;

%               d_period= 2*pi/x(4) % the result should be consitent with the value lamda/(2n)
           end
           
 %%%%
 % Print the images
 cd (folder);
mkdir('.\_FringesResults'); 
FigureName=fullfile ('.\_FringesResults\','_PeaksCentre');% get full path
figure (6)
fig = gcf;
fig.Color = 'white'; % set the background figure color
fig.InvertHardcopy = 'off';
iptsetpref('ImshowBorder','tight'); % Figures without any borders 
print(fig,FigureName,'-dtiff','-r100');

figure (7)
FigureName=fullfile ('.\_FringesResults\','_CentreHeight-time');% get full path
fig = gcf;
fig.Color = 'white'; % set the background figure color
fig.InvertHardcopy = 'off';
iptsetpref('ImshowBorder','tight'); % Figures without any borders 
print(fig,FigureName,'-dtiff','-r100');
close all;
%% this the section to plot the cross-section of the droplet and reconstruct the drop profile 
%!! From the refence of the contact line
cd (oldFolder);
 
smooth_points_r=6; % for tracing cross section
smooth_points_r_head=2;
smooth_points_r_tail=2;
threshold_peaks_r=4;  %set the criteria to distinguish peaks and valleys

 backwards_index=4; % use it as the final peak/valley 'backwards_index=0' is not stable when tracing the edge of the bank
 index_analysis_frame0=length(h_c)-backwards_index;% h_c only associates with bright and dark ring at the centre
 R_c_profile=zeros(index_analysis_frame0,1);
 h_c_profile=zeros(index_analysis_frame0,1);
 CA_profile=zeros(index_analysis_frame0,1);
 drying_time_fraction=zeros(index_analysis_frame0,1);
 profileCell_r={}; % container for profiles of different frames, for the radius
 profileCell_h={}; % container for profiles of different frames, for the height
 %frame_from_dried=round(h_c(index_analysis_frame,1)); % reversed frames which corresonds to the index_analysis_frame above, skip frame considered
%  frame_from_start_real=analysis_start-(frame_from_dried-1)*(skip_frames+1);
 % drying_time_fraction=round(100*(frame_from_start_real-drying_start_frame)/drying_time_in_frame)/100;
 
 %for drop_frame=analysis_start:-(skip_frames+1):analysis_end
 count=1;
 
 % For plot figures
colormark=[0    0.4470    0.7410;
    0.8500    0.3250    0.0980;
    0.9290    0.6940    0.1250;
    0.4940    0.1840    0.5560;
    0.4660    0.6740    0.1880;
    0.3010    0.7450    0.9330];
colorCount=0;
 for t0= index_analysis_frame0:-1:2 
 
 h0=h_c(t0,2);% centre info
 frame_from_dried=round(h_c(t0,1)); % reversed frames which corresonds to the index_analysis_frame above, skip frame considered
 frame_from_start_real=analysis_start-(frame_from_dried-1)*(skip_frames+1);
 drying_time_fraction(t0)=round(1000*(frame_from_start_real-drying_start_frame)/drying_time_in_frame)/1000;
 k=frame_from_dried;
 I=I_3D(:,:,frame_from_dried); % this corresponds to infor of the frame_from_start_real
  
    figure(9)
    cla
    imagesc(I);
    colormap gray; 
    axis equal tight;
    hold on    
            
    th=0:0.01:2*pi; %%%%Angle
    xy_c=[xc_c2(k)+R_c2(k)*cos(th) ; yc_c2(k)+R_c2(k)*sin(th)]';

    plot(xy_c(:,1),xy_c(:,2),'-r');
    plot(xc_c2(k),yc_c2(k),'xr');
    
      %%% Horizontal Lines across the centre
      theta=0; % 0.5*pi for vertical lines
       % the end points to define the line
      x=[xc_c2(k)-R_c2(k)*cos(theta) xc_c2(k)+R_c2(k)*cos(theta)];
      y=[yc_c2(k)-R_c2(k)*sin(theta) yc_c2(k)+R_c2(k)*sin(theta)];

          
        line(x,y,'Color','y'); % a line that cross the centre and the input point and within the CL
   
        [cx cy intensity]=improfile(I,x,y);
        % improfile retrieves the intensity values of pixels along a line 
        % interpolation 'nearest', I have checked that intensity(i)=I(round(cy(i)),round(cx(i)))
        r=((cx-x(1)).^2+(cy-y(1)).^2).^(1/2);
        length_c=((x(2)-x(1)).^2+(y(2)-y(1)).^2).^(1/2)/2; % in case of circle, equals R_c
        r_c=r-length_c;% shift drop centre to zero x-coordinate
        
        r0=0;%((xc_c-x(1)).^2+(yc_c-y(1)).^2).^(1/2)-length_c;% centre, r0 should be close to zero
        I0=I(round(yc_c2(k)),round(xc_c2(k)));
%         h0; % um centre
%         frame_from_dried;
        
        cross_r=[cx cy r_c intensity];
      %  if ~mod(length(temp),2) 
        cross_r=[cross_r; xc_c2(k) yc_c2(k) r0 I0]; % add [r0 I0] so that centre is well defined
     %   end
        cross_r=sortrows(cross_r,3);
         
        figure(10)
        cla
        hold on
        plot(cross_r(:,3), cross_r(:,4),'b-o');
         xlabel('r /pixels');
         ylabel('Intensity');
         axis tight
         
  %Intial plot of smooth; To find the devideFraction
        
        y1=smooth(double(cross_r(:,4)),smooth_points_r,'lowess'); % the second parameter is important. Adjust it according to your image contrast and noise level
        y2=smooth(double(cross_r(:,4)),0.5,'lowess');
        intensity=y1-y2;
        plot(cross_r(:,3),intensity,'r');
       
        [maxtab, mintab] = peakdet(intensity, threshold_peaks_r, cross_r(:,3)); % keep the second parameter as low as possible. but be aware of noises. 
        index_head_divide=find(cross_r(:,3)==mintab(2,1)); % find the 2 valley index in cross_r, this point as the devide point for smooth_points_r_head
        index_tail_divide=find(cross_r(:,3)==maxtab(end-2,1)); % find the 3-to-end peak index in cross_r, this point as the devide point for smooth_points_r_head
        hold on
         plot(mintab(:,1), mintab(:,2), 'g*');
         plot(maxtab(:,1), maxtab(:,2), 'r*');  
%          plot(cross_r(index_0,3),intensity(index_0),'ro');
        % Check the figure to decide the left and right divide points
        NumCross_r=length (cross_r);
        if length (mintab)>3 && length (maxtab)>3
        left_divide_point=3;
        right_divide_point=3;
        index_head_divide=find(cross_r(:,3)==mintab(left_divide_point,1)); % find the 2 valley index in cross_r, this point as the devide point for smooth_points_r_head
        index_tail_divide=find(cross_r(:,3)==maxtab(end-right_divide_point-1,1)); % find the 3-to-end peak index in cross_r, this point as the devide point for smooth_points_r_head
        else
         DevideFraction=0.2; % this fraction needs to be checked
         index_head_divide=round(DevideFraction*NumCross_r);
         index_tail_divide=NumCross_r-round(DevideFraction*NumCross_r)+1;
        end
 %%% devide the cross_r into two parts: the central part and the edge part 
 %%% This is for separate smooth when the edge part has very dense fringes 
 %%% Comment this if not needed
         if index_head_divide<round(0.3*NumCross_r)&& (NumCross_r-index_tail_divide)<round(0.3*NumCross_r)
             HeadCross_r=cross_r(1:index_head_divide,:);
             TailCross_r=cross_r(index_tail_divide:end,:);
             MiddleCross_r=cross_r(index_head_divide+1:index_tail_divide-1,:);

             HeadSmooth=smooth(double(HeadCross_r(:,4)),smooth_points_r_head,'lowess'); % the second parameter is important. Adjust it according to your image contrast and noise level
             TailSmooth=smooth(double(TailCross_r(:,4)),smooth_points_r_tail,'lowess'); % the second parameter is important. Adjust it according to your image contrast and noise level
             MiddleSmooth=smooth(double(MiddleCross_r(:,4)),smooth_points_r,'lowess'); % the second parameter is important. Adjust it according to your image contrast and noise level
             y1=[HeadSmooth;MiddleSmooth;TailSmooth];
         end
%%%        y1=smooth(double(cross_r(:,4)),smooth_points_r,'lowess'); % the second parameter is important. Adjust it according to your image contrast and noise level

         figure (13)
         cla
         plot(cross_r(:,3),y1,'r');
         y2=smooth(y1,0.5,'lowess');
%          figure (14)
%          cla
%          plot(cross_r(:,3),y2,'r');      

%         figure (11)
%         cla
%         y1=smooth(double(cross_r(:,4)),smooth_points_r,'lowess'); % the second parameter is important. Adjust it according to your image contrast and noise level
%         y2=smooth(double(cross_r(:,4)),0.5,'lowess');
%         figure (12)
%         cla
%         plot(cross_r(:,3),y2,'r');
        intensity=y1-y2;
%         intensity=y1;
        figure (15)
        cla
        plot(cross_r(:,3),intensity,'r');
        hold on
        
        height_r=[];
        profile=[];
        index_0=find(abs(cross_r(:,3))==0);
        if t0/index_analysis_frame0<=0.3
          threshold_peaks_r=8;    
        end
          
                 [maxtab, mintab] = peakdet(intensity, threshold_peaks_r, cross_r(:,3)); % keep the second parameter as low as possible. but be aware of noises. 
                 % delete the over-counted centre values afterwards
                 if ~isempty(mintab) & ~isempty(maxtab)
                     fringes_r=sortrows([maxtab;mintab;[cross_r(index_0,3) intensity(index_0)]],1);
             
                     plot(mintab(:,1), mintab(:,2), 'g*');
                     plot(maxtab(:,1), maxtab(:,2), 'r*');  
                     plot(cross_r(index_0,3),intensity(index_0),'ro');
                     
                     xlabel('r /pixels');
                     ylabel('Intensity');
                     axis tight      
                      
                     region_beyond=10; % the region beyond the centre fringe, check this value        
                     index_ROI=find(fringes_r(:,1)<-region_beyond);  
                       
                        for i=1:length(index_ROI)
                          height_r=[h0-i*lambda/(4*n_index);height_r];
                        end
                        profile=[[fringes_r(index_ROI,1) height_r];[r0 h0]]; %[r h]
                                       
                        index_ROI=find(fringes_r(:,1)>region_beyond); 
                        for i=index_ROI(1):index_ROI(end)
                          profile=[profile;[fringes_r(i) -(i-index_ROI(1)+1)*lambda/(4*n_index)+h0]];
                        end              
                      
       
                        % fit the profile with circle
                        circ_fit=[profile ones(length(profile),1)]\[-(profile(:,1).^2+profile(:,2).^2)];
                        xc_c_profile = -.5*circ_fit(1);
                        yc_c_profile = -.5*circ_fit(2);
                        R_c_profile (t0)  =  sqrt((circ_fit(1)^2+circ_fit(2)^2)/4-circ_fit(3));
                        h_c_profile(t0)=R_c_profile(t0)+yc_c_profile;
                        CA_profile(t0)=2*atan(h_c_profile(t0)/R_c_profile (t0))*(180/pi); % contact angle tan(2CA)=h/R
                  
                        th_profile=pi/4:0.001:pi/2+pi/4; %%%%Angle              
                        xy_c_profile=[xc_c_profile+R_c_profile(t0)*cos(th_profile) ; yc_c_profile+R_c_profile(t0)*sin(th_profile)]';
                        id=find((xy_c_profile(:,1)>=-length_c-2)&(xy_c_profile(:,1)<=length_c+2)&(xy_c_profile(:,2)>0));
                        figure(11)
                        if t0== index_analysis_frame0
                            cla
                        end
                        if mod(t0,2)==0
                        hold on
                        colorCount=colorCount+1;
                        plot(profile(:,1),profile(:,2),'*','MarkerEdgeColor', 'k');                  
                        plot(xy_c_profile(id,1),xy_c_profile(id,2),'-','color',[0,0,0]+0.5,'LineWidth',1.2)
                        
                        xlabel('r /pixels');
                        ylabel('Height/\mum');
                        axis tight 
                        end
                        
                 end
                 profileCell_r{t0}=profile(:,1);
                 profileCell_h{t0}=profile(:,2);
%                   pause;
end
pause;
figure (11)
legend('off')
fig = gcf;
fig.Color = 'white'; % set the background figure color
fig.InvertHardcopy = 'off';
iptsetpref('ImshowBorder','tight'); % Figures without any borders 
cd (folder)
FigureName=fullfile ('.\_FringesResults\','_Profile_R-h');
print(fig,FigureName,'-dtiff','-r100' );

% Plot the contact angle with time
figure (1)
cla
plot(drying_time_fraction(2:end),CA_profile (2:end),'-*')
% set(gca, 'XLim', [1, 2000],'YLim', [0.1, 2])
xlabel('t/t_{dry}');
ylabel('CA / degree');

fig = gcf;
fig.Color = 'white'; % set the background figure color
fig.InvertHardcopy = 'off';
iptsetpref('ImshowBorder','tight'); % Figures without any borders 
FigureName=fullfile ('.\_FringesResults\','_CA-t');
print(fig,FigureName,'-dtiff','-r100' );

%% this the section to plot the cross-section of the droplet and reconstruct the drop profile 
%!! From the refence of the centre point
 cd (oldFolder);
 
smooth_points_r=6; % for tracing cross section
smooth_points_r_head=2;
smooth_points_r_tail=2;
threshold_peaks_r=4;  %set the criteria to distinguish peaks and valleys

 backwards_index=4; % use it as the final peak/valley 'backwards_index=0' is not stable when tracing the edge of the bank
 index_analysis_frame0=length(h_c)-backwards_index;% h_c only associates with bright and dark ring at the centre
 R_c_profile=zeros(index_analysis_frame0,1);
 h_c_profile=zeros(index_analysis_frame0,1);
 CA_profile=zeros(index_analysis_frame0,1);
 drying_time_fraction=zeros(index_analysis_frame0,1);
 profileCell_r={}; % container for profiles of different frames, for the radius
 profileCell_h={}; % container for profiles of different frames, for the height
 %frame_from_dried=round(h_c(index_analysis_frame,1)); % reversed frames which corresonds to the index_analysis_frame above, skip frame considered
%  frame_from_start_real=analysis_start-(frame_from_dried-1)*(skip_frames+1);
 % drying_time_fraction=round(100*(frame_from_start_real-drying_start_frame)/drying_time_in_frame)/100;
 
 %for drop_frame=analysis_start:-(skip_frames+1):analysis_end
 count=1;
 
 % For plot figures
colormark=[0    0.4470    0.7410;
    0.8500    0.3250    0.0980;
    0.9290    0.6940    0.1250;
    0.4940    0.1840    0.5560;
    0.4660    0.6740    0.1880;
    0.3010    0.7450    0.9330];
colorCount=0;
 for t0= index_analysis_frame0:-1:2 
 
 h0=h_c(t0,2);% centre info
 frame_from_dried=round(h_c(t0,1)); % reversed frames which corresonds to the index_analysis_frame above, skip frame considered
 frame_from_start_real=analysis_start-(frame_from_dried-1)*(skip_frames+1);
 drying_time_fraction(t0)=round(1000*(frame_from_start_real-drying_start_frame)/drying_time_in_frame)/1000;
 k=frame_from_dried;
 I=I_3D(:,:,frame_from_dried); % this corresponds to infor of the frame_from_start_real
  
    figure(9)
    cla
    imagesc(I);
    colormap gray; 
    axis equal tight;
    hold on    
            
    th=0:0.01:2*pi; %%%%Angle
    xy_c=[xc_c2(k)+R_c2(k)*cos(th) ; yc_c2(k)+R_c2(k)*sin(th)]';

    plot(xy_c(:,1),xy_c(:,2),'-r');
    plot(xc_c2(k),yc_c2(k),'xr');
    
      %%% Horizontal Lines across the centre
      theta=0; % 0.5*pi for vertical lines
       % the end points to define the line
      x=[xc_c2(k)-R_c2(k)*cos(theta) xc_c2(k)+R_c2(k)*cos(theta)];
      y=[yc_c2(k)-R_c2(k)*sin(theta) yc_c2(k)+R_c2(k)*sin(theta)];

          
        line(x,y,'Color','y'); % a line that cross the centre and the input point and within the CL
   
        [cx cy intensity]=improfile(I,x,y);
        % improfile retrieves the intensity values of pixels along a line 
        % interpolation 'nearest', I have checked that intensity(i)=I(round(cy(i)),round(cx(i)))
        r=((cx-x(1)).^2+(cy-y(1)).^2).^(1/2);
        length_c=((x(2)-x(1)).^2+(y(2)-y(1)).^2).^(1/2)/2; % in case of circle, equals R_c
        r_c=r-length_c;% shift drop centre to zero x-coordinate
        
        r0=0;%((xc_c-x(1)).^2+(yc_c-y(1)).^2).^(1/2)-length_c;% centre, r0 should be close to zero
        I0=I(round(yc_c2(k)),round(xc_c2(k)));
%         h0; % um centre
%         frame_from_dried;
        
        cross_r=[cx cy r_c intensity];
      %  if ~mod(length(temp),2) 
        cross_r=[cross_r; xc_c2(k) yc_c2(k) r0 I0]; % add [r0 I0] so that centre is well defined
     %   end
        cross_r=sortrows(cross_r,3);
         
        figure(10)
        cla
        hold on
        plot(cross_r(:,3), cross_r(:,4),'b-o');
         xlabel('r /pixels');
         ylabel('Intensity');
         axis tight
         
  %Intial plot of smooth; To find the devideFraction
        
        y1=smooth(double(cross_r(:,4)),smooth_points_r,'lowess'); % the second parameter is important. Adjust it according to your image contrast and noise level
        y2=smooth(double(cross_r(:,4)),0.5,'lowess');
        intensity=y1-y2;
        plot(cross_r(:,3),intensity,'r');
       
        [maxtab, mintab] = peakdet(intensity, threshold_peaks_r, cross_r(:,3)); % keep the second parameter as low as possible. but be aware of noises. 
        index_head_divide=find(cross_r(:,3)==mintab(2,1)); % find the 2 valley index in cross_r, this point as the devide point for smooth_points_r_head
        index_tail_divide=find(cross_r(:,3)==maxtab(end-2,1)); % find the 3-to-end peak index in cross_r, this point as the devide point for smooth_points_r_head
        hold on
         plot(mintab(:,1), mintab(:,2), 'g*');
         plot(maxtab(:,1), maxtab(:,2), 'r*');  
%          plot(cross_r(index_0,3),intensity(index_0),'ro');
        % Check the figure to decide the left and right divide points
        NumCross_r=length (cross_r);
        if length (mintab)>3 && length (maxtab)>3
        left_divide_point=3;
        right_divide_point=3;
        index_head_divide=find(cross_r(:,3)==mintab(left_divide_point,1)); % find the 2 valley index in cross_r, this point as the devide point for smooth_points_r_head
        index_tail_divide=find(cross_r(:,3)==maxtab(end-right_divide_point-1,1)); % find the 3-to-end peak index in cross_r, this point as the devide point for smooth_points_r_head
        else
         DevideFraction=0.2; % this fraction needs to be checked
         index_head_divide=round(DevideFraction*NumCross_r);
         index_tail_divide=NumCross_r-round(DevideFraction*NumCross_r)+1;
        end
 %%% devide the cross_r into two parts: the central part and the edge part 
 %%% This is for separate smooth when the edge part has very dense fringes 
 %%% Comment this if not needed
         if index_head_divide<round(0.3*NumCross_r)&& (NumCross_r-index_tail_divide)<round(0.3*NumCross_r)
             HeadCross_r=cross_r(1:index_head_divide,:);
             TailCross_r=cross_r(index_tail_divide:end,:);
             MiddleCross_r=cross_r(index_head_divide+1:index_tail_divide-1,:);

             HeadSmooth=smooth(double(HeadCross_r(:,4)),smooth_points_r_head,'lowess'); % the second parameter is important. Adjust it according to your image contrast and noise level
             TailSmooth=smooth(double(TailCross_r(:,4)),smooth_points_r_tail,'lowess'); % the second parameter is important. Adjust it according to your image contrast and noise level
             MiddleSmooth=smooth(double(MiddleCross_r(:,4)),smooth_points_r,'lowess'); % the second parameter is important. Adjust it according to your image contrast and noise level
             y1=[HeadSmooth;MiddleSmooth;TailSmooth];
         end
%%%        y1=smooth(double(cross_r(:,4)),smooth_points_r,'lowess'); % the second parameter is important. Adjust it according to your image contrast and noise level

         figure (13)
         cla
         plot(cross_r(:,3),y1,'r');
         y2=smooth(y1,0.5,'lowess');
%          figure (14)
%          cla
%          plot(cross_r(:,3),y2,'r');      

%         figure (11)
%         cla
%         y1=smooth(double(cross_r(:,4)),smooth_points_r,'lowess'); % the second parameter is important. Adjust it according to your image contrast and noise level
%         y2=smooth(double(cross_r(:,4)),0.5,'lowess');
%         figure (12)
%         cla
%         plot(cross_r(:,3),y2,'r');
        intensity=y1-y2;
%         intensity=y1;
        figure (15)
        cla
        plot(cross_r(:,3),intensity,'r');
        hold on
        
        height_r=[];
        profile=[];
        index_0=find(abs(cross_r(:,3))==0);
        if t0/index_analysis_frame0<=0.3
          threshold_peaks_r=8;    
        end
          
                 [maxtab, mintab] = peakdet(intensity, threshold_peaks_r, cross_r(:,3)); % keep the second parameter as low as possible. but be aware of noises. 
                 % delete the over-counted centre values afterwards
                 if ~isempty(mintab) & ~isempty(maxtab)
                     fringes_r=sortrows([maxtab;mintab;[cross_r(index_0,3) intensity(index_0)]],1);
             
                     plot(mintab(:,1), mintab(:,2), 'g*');
                     plot(maxtab(:,1), maxtab(:,2), 'r*');  
                     plot(cross_r(index_0,3),intensity(index_0),'ro');
                     
                     xlabel('r /pixels');
                     ylabel('Intensity');
                     axis tight      
                      
                     region_beyond=10; % the region beyond the centre fringe, check this value        
                     index_ROI=find(fringes_r(:,1)<-region_beyond);  
                       
                        for i=1:length(index_ROI)
                          height_r=[h0-i*lambda/(4*n_index);height_r];
                        end
                        profile=[[fringes_r(index_ROI,1) height_r];[r0 h0]]; %[r h]
                                       
                        index_ROI=find(fringes_r(:,1)>region_beyond); 
                        for i=index_ROI(1):index_ROI(end)
                          profile=[profile;[fringes_r(i) -(i-index_ROI(1)+1)*lambda/(4*n_index)+h0]];
                        end              
                      
       
                        % fit the profile with circle
                        circ_fit=[profile ones(length(profile),1)]\[-(profile(:,1).^2+profile(:,2).^2)];
                        xc_c_profile = -.5*circ_fit(1);
                        yc_c_profile = -.5*circ_fit(2);
                        R_c_profile (t0)  =  sqrt((circ_fit(1)^2+circ_fit(2)^2)/4-circ_fit(3));
                        h_c_profile(t0)=R_c_profile(t0)+yc_c_profile;
                        CA_profile(t0)=2*atan(h_c_profile(t0)/R_c_profile (t0))*(180/pi); % contact angle tan(2CA)=h/R
                  
                        th_profile=pi/4:0.001:pi/2+pi/4; %%%%Angle              
                        xy_c_profile=[xc_c_profile+R_c_profile(t0)*cos(th_profile) ; yc_c_profile+R_c_profile(t0)*sin(th_profile)]';
                        id=find((xy_c_profile(:,1)>=-length_c-2)&(xy_c_profile(:,1)<=length_c+2)&(xy_c_profile(:,2)>0));
                        figure(11)
                        if t0== index_analysis_frame0
                            cla
                        end
                        if mod(t0,2)==0
                        hold on
                        colorCount=colorCount+1;
                        plot(profile(:,1),profile(:,2),'*','MarkerEdgeColor', 'k');                  
                        plot(xy_c_profile(id,1),xy_c_profile(id,2),'-','color',[0,0,0]+0.5,'LineWidth',1.2)
                        
                        xlabel('r /pixels');
                        ylabel('Height/\mum');
                        axis tight 
                        end
                        
                 end
                 profileCell_r{t0}=profile(:,1);
                 profileCell_h{t0}=profile(:,2);
%                   pause;
end
pause;
figure (11)
legend('off')
fig = gcf;
fig.Color = 'white'; % set the background figure color
fig.InvertHardcopy = 'off';
iptsetpref('ImshowBorder','tight'); % Figures without any borders 
cd (folder)
FigureName=fullfile ('.\_FringesResults\','_Profile_R-h');
print(fig,FigureName,'-dtiff','-r100' );

% Plot the contact angle with time
figure (1)
cla
plot(drying_time_fraction(2:end),CA_profile (2:end),'-*')
% set(gca, 'XLim', [1, 2000],'YLim', [0.1, 2])
xlabel('t/t_{dry}');
ylabel('CA / degree');

fig = gcf;
fig.Color = 'white'; % set the background figure color
fig.InvertHardcopy = 'off';
iptsetpref('ImshowBorder','tight'); % Figures without any borders 
FigureName=fullfile ('.\_FringesResults\','_CA-t');
print(fig,FigureName,'-dtiff','-r100' );


  
%       %%
%        y=zeros(pixels_traced,1);
%        frame_from_dried=800; % reversed frames
%        frame_from_start_real=analysis_start-(frame_from_dried-1)*(skip_frames+1)
%        drying_time_fraction=round(100*(frame_from_start_real-drying_start_frame)/drying_time_in_frame)/100
%        
%        y(:)=H_r_t(frame_from_dried,2,:);
%        figure(11)       
%        plot(H_t0(:,3),y,'o')
% %% save file
%   cd(folder)
%   
%     save I_3D.mat I_3D
%     save H_r_t.mat H_r_t
%     save H_t0.mat H_t0
%     cd ..
       
       