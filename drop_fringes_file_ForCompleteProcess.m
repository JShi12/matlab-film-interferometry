%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% files
% file example: fringes
% Beth's OLED film drying
% Example IOP\ sample: anisole/methyl benezoate 
fileformat ='png';
folder='H:\M-Matlab codes\Film interferometry\Surfactants\Images';
filetitle='Frame';

background_file_start=1; % with flat surface before drop landing 
background_file_end=1080;% the frame where inner film disappears
deposit_file=1080; % the frame where inner film disappears

drying_start_frame=1; % whole drop landed
drying_end_frame=1080; % end of drying


inverse_frame_start=1079;% The frame where there is a very small inner ring
inverse_frame_end=2;% end frame of your interest
frame_RingForm=686;% frame when a ring forms
skip_frames=0; % take frame+skip_frames+1, no skip: 0



%%% resolutions
resolution=0.40;%um/pixel
frame_rate=1000/(skip_frames+1);%fps 125 for test_01; water case 100fps, 50ethanol/water: \test_image_1: 3000fps. it has to be carefully done as camera could skip save frames
frame_dt=1/frame_rate;%time between neighbor frames /s

t_dry=(drying_end_frame-drying_start_frame)*frame_dt*1000; % total drying time in ms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% deposit shapes
circ_dot=1;

%optics
lambda=0.465; %um blue LED thorlab 470nm
n_index=1.561; % anisole perhaps