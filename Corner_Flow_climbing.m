     
clc
clear
close all

%% Data configuration -- need to change
listing = dir('/Volumes/Judy_SoilC_2/Bacteria_climbing/Climbing_with_bacteria/34oC_incubator/2chambers_30/June25_WT_OD0p23_0p004NBDG/TimeLapse_9030_55ul_16hr_1min');
nb_all_img = 995; %numel(listing)-2;
dt = 61/60; % minutes 
t0 = 0; % minutes
Intensity_crit = 0.05; %% ?????? has to tune this value


% Width_chamber = 620; % ???? in pixels need to confirm this number from images 
% magnification =1/2; 

pixel_size = 16/sqrt(3)/836; %mm
%3.45/10^3/magnification; % depend on camera magnification ratio, 1:1 ratio is 3.45um per pixel 

output_video_name = 'Track_boundaryright.avi'; %s??????
%% configuration 

gap = 10;
gap_plot = 1;
T_all = [t0:dt*gap:t0+dt*(nb_all_img-1)];
T_all_nogap = [t0:dt:t0+dt*(nb_all_img-1)];

%% save the images with imposed boundary to a video 
workingDir = 'processed2';
mkdir(workingDir)

outputVideo = VideoWriter(fullfile(workingDir,output_video_name));
outputVideo.FrameRate = 10;  %% output Framerate , can be changed ??????
open(outputVideo)

%% crop images
% waitfor(msgbox('Crop the region of interest'));
A = imread([listing(nb_all_img+2).folder,'/',listing(nb_all_img+2).name]);%nb_all_img
A0 = A;
imshow(A)

A_gray = rgb2gray(A);
A_intensity = A_gray ;
A_intensity = medfilt2(A_intensity);

figure
imshow(A_intensity)

% rect = getrect % [xmin ymin width height]

rect = [ 208          50         915        1464 ]; % right chamber
A = imcrop(A,rect);

%% Crop the region for background intensity correction 
% waitfor(msgbox('Crop the region for background'));
% A = imread(fname, 1);
imshow(A0)
% rect_bkg = getrect % [xmin ymin width height]
rect_bkg = [208          50         915        1464]; %right chamber 

A_bkg0 = imcrop(A0,rect_bkg);
A_gray_bkg0 = rgb2gray(A_bkg0);
BKG_intensity_0 = mean(mean(A_gray_bkg0))/255;
BKG_intensity_0 = 0

%% Crop the region for tracking the left boundary of the bacteria 
% waitfor(msgbox('Crop the region for tracking the left side of the bacteria (polygon)'));
% A = imread(fname, 1);
imshow(A)
% [rect_left,xi2,yi2] = roipoly; % [xmin ymin width height]

% xi2 = [109   145    38    23    38    88   109]'; %right chamber
% yi2 = [1065        1158        1299         988          59          56        1065]';

xi2 = [  14     8   115   133    14]';
yi2 = [1161         194         206        1146        1161]';
%% Crop the region for tracking the right boundary of the bacteria 
% waitfor(msgbox('Crop the region for tracking the right side  of the bacteria (polygon)'));
% A = imread(fname, 1);
imshow(A)
[rect_right,xi3,yi3] = roipoly; % [xmin ymin width height]
% 
% xi3 = [798   667   730   748   822   798]'; %right chamber
% yi3 = [1203        1092         976          17          11        1203]';

xi3 = [ 583   703   897   915   583 ]';
yi3 = [1254          14          17        1275        1254]';
% %% Crop the region for tracking the middle region of the bacteria 
% % waitfor(msgbox('Crop the region for tracking the middle edge of the bacteria (polygon)'));
% % A = imread(fname, 1);
% imshow(A)
% % [rect_mid,xi4,yi4] = roipoly; % [xmin ymin width height]
% 
% xi4 = [ 316   337   517   502   316]';
% yi4 = [ 1074  17  23  1083  1074]';

%             figure
%             imshow(A)
%             hold on
% 
%             plot(xi2,yi2, 'r', 'LineWidth', 2)

%         A_cut = imcrop(A,rect);
% A_cut = A.* uint8(rect);

%% find the left and right boundary of the chamber

leftside_x = [44 5]; 
leftside_y = [1460 8];

rightside_x = [909 876]; 
rightside_y = [1443 1];

%% image process 

bnd_pic_ind = [];
bnd_pic_ind_up = [];
below_cloud_int = [];
Left_jet_index = [];
Right_jet_index = [];
Mid_jet_index = [];
for k = 1:gap:nb_all_img
    
    %%%  Step 1: read iamges %%%%%
    A0 = imread([listing(k+2).folder,'/',listing(k+2).name]);
    A = imcrop(A0,rect);
    A_bkg = imcrop(A0,rect_bkg);


    %%%  Step 2: change to gray scale %%%%%
    A_gray = rgb2gray(A);
    A_gray_bkg = rgb2gray(A_bkg);
    BKG_intensity(k) = mean(mean(A_gray_bkg))/255;
    A_intensity = A_gray ;
    
    A_intensity = medfilt2(A_intensity);
    
    %%% calculate the intensity 
    
    I_tot(k) = sum(sum(A_intensity));
    vertical(:,k) = mean(A_intensity,2);
    horizontal(:,k) = mean(A_intensity,1);

%     figure
%     imshow(A_intensity) 

    
        %%  Step 3: make image to binary images %%%%%
            BW = imbinarize(A_intensity,Intensity_crit+BKG_intensity(k) -BKG_intensity_0);
            
            BW = bwareaopen(BW,2000); % remove objects less than 50 pixels
%             figure;
%             imshow(BW)


        
         %% Step 4: find the boundary of the cloud (the largest object)
         
            [Bdary,Bdary_label,nb_objects, Bdary_A] = bwboundaries(BW,'noholes');

            Bdary_size = cellfun(@numel, Bdary);
            [Bdary_dots, Bdary_ID] = max(Bdary_size);
            [B,I_bdary_size] = sort(Bdary_size );
            
            %% find the largest region as the bacteria region   
            Bdary_bac_x = Bdary{Bdary_ID}(:,2);
            Bdary_bac_y = Bdary{Bdary_ID}(:,1);
            Bac_BW = poly2mask(Bdary_bac_x,Bdary_bac_y,rect(4)+1,rect(3)+1);

            %% calulate average bacteria intensity 
             bac_tot(k) = sum(A_intensity(find(Bac_BW>0)));   
             area_bac(k) = numel(find(Bac_BW>0));
             I_ave(k) = bac_tot(k)/area_bac(k);
        %% Plot the boundary of the clound in the original image
        
        h = figure(3);  
%             set(h,'visible','off')
        imshow(A); hold on
        line(leftside_x,leftside_y,'color','c','linewidth',2)
        line(rightside_x,rightside_y,'color','c','linewidth',2)
%             
        if rem(k,gap_plot)==0

%             plot(xi2,yi2, 'Color', [17 17 17]/255, 'LineWidth', 2) 
            hold on
        end
            
        if isempty(Bdary_ID)==0
            if rem(k,gap_plot)==0
            plot(Bdary{Bdary_ID}(:,2), Bdary{Bdary_ID}(:,1), 'k', 'LineWidth', 2);
            end
            

%             for ii2 = 1:numel(Bdary)
%                  Bdary_bac_x = cat(1,Bdary_bac_x,Bdary{ii2}(:,2));
%                  Bdary_bac_y = cat(1,Bdary_bac_y,Bdary{ii2}(:,1));
%                  if rem(k,gap_plot)==0
%                  plot(Bdary{ii2}(:,2), Bdary{ii2}(:,1), 'k', 'LineWidth', 2);
%                  end
%             end
            
%             %% find the boundary inside the left, right, and middle region
%             
%             left_jet_ind = inpolygon(Bdary_bac_x, Bdary_bac_y,xi2,yi2);
%             right_jet_ind = inpolygon(Bdary_bac_x, Bdary_bac_y,xi3,yi3);
%             mid_jet_ind = inpolygon(Bdary_bac_x, Bdary_bac_y,xi4,yi4);
            
            %% find the boundary inside the left, right, and middle region

            left_jet_ind = inpolygon(Bdary_bac_x, Bdary_bac_y,xi2,yi2);
            right_jet_ind = inpolygon(Bdary_bac_x, Bdary_bac_y,xi3,yi3);
            
            %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if sum(left_jet_ind )>0
%                 if rem(k,gap_plot)==0
%                     plot(Bdary_bac_x(left_jet_ind) , Bdary_bac_y(left_jet_ind) , 'r.', 'LineWidth', 2);
%                 end
                    left_jet_x{k} = Bdary_bac_x(left_jet_ind);
                    left_jet_y{k} = Bdary_bac_y(left_jet_ind);
                    ymax_LeftJet(k) = min(left_jet_y{k});
                    cloud_y_sort = sort(left_jet_y{k});
                    %% calculate the width of the left jet 
                    width_Ljet = [];
                    for y_ii = 1:numel(cloud_y_sort)
                        index_jet = find(left_jet_y{k} ==cloud_y_sort(y_ii));
                        wall_left{k}(y_ii) = interp1(leftside_y,leftside_x,cloud_y_sort(y_ii));
                        bac_left{k}(y_ii) = max(left_jet_x{k}(index_jet));
                        width_Ljet{k}(y_ii) = bac_left{k}(y_ii)-wall_left{k}(y_ii);
%                         if(numel(index_jet)==2)
%                             width_Ljet =  cat(1,width_Ljet,max(left_jet_x{k}(index_jet))-min(left_jet_x{k}(index_jet)));                   
%                         end
                    end
                    cloud_left{k} = cloud_y_sort;
                    if rem(k,gap_plot)==0
                        plot(wall_left{k}, cloud_left{k} , 'r.', 'LineWidth', 2);
                        plot(bac_left{k}, cloud_left{k} , 'r.', 'LineWidth', 2);
                    end
                
                        width_Ljet_mean(k) = mean(width_Ljet{k});     
                        width_Ljet_SD(k) = std(width_Ljet{k}); 
                        
                        Left_jet_index  = cat(1, Left_jet_index  ,k);

            end
            if sum(right_jet_ind )>0
%                 if rem(k,gap_plot)==0
%                     plot(Bdary_bac_x(right_jet_ind) , Bdary_bac_y(right_jet_ind) , 'b.', 'LineWidth', 2);
%                 end
                    right_jet_x{k} = Bdary_bac_x(right_jet_ind);
                    right_jet_y{k} = Bdary_bac_y(right_jet_ind);
                    ymax_RightJet(k) = min(right_jet_y{k});
                    cloud_y_sort = sort(right_jet_y{k});
                    %% calculate the width of the right jet 
                    width_Rjet = [];
                    for y_ii = 1:numel(cloud_y_sort)
                        index_jet = find(right_jet_y{k} ==cloud_y_sort(y_ii)); 
                        wall_right{k}(y_ii) = interp1(rightside_y,rightside_x,cloud_y_sort(y_ii));
                        bac_right{k}(y_ii) = min(right_jet_x{k}(index_jet));
                        width_Rjet{k}(y_ii) = wall_right{k}(y_ii)-bac_right{k}(y_ii);

                    end
                        width_Rjet_mean(k) = mean(width_Rjet{k});     
                        width_Rjet_SD(k) = std(width_Rjet{k}); 
                        
                    cloud_right{k} = cloud_y_sort;
                    if rem(k,gap_plot)==0
                        plot(wall_right{k}, cloud_right{k} , 'b.', 'LineWidth', 2);                    
                        plot(bac_right{k}, cloud_right{k} , 'b.', 'LineWidth', 2);
                    end
                        
                        Right_jet_index  = cat(1, Right_jet_index  ,k);
            end
            
            
            
            
             
        end
            

%             frame = getframe(gcf);
%             writeVideo(outputVideo,frame)
            
               set(gca,'xtick',[],'ytick',[])
               set(gca,'LooseInset',get(gca,'TightInset'));
    %             set(gcf, 'PaperSize', [4 2]);
    %             fig = gcf
    %             set(fig,'PaperPositionMode','auto')

                saveas(gcf,...
                    [pwd,'/WT_nogfp_30/','t_',num2str(k),'.tif']);
                 


            end
            

%     set(h_pdf,'PaperPositionMode','Manual')
%     set(h_pdf,'PaperUnits','inches')
%     set(h_pdf,'PaperSize',[8 6])
%     set(h_pdf,'PaperPosition',[0 0 8 6])
%     saveas(h_pdf, ['Img_intensity_max'], 'pdf')
%     
    close(outputVideo)
 
     %% plot the height of the cloud as a function of time 
    
    figure; hold on
    set(gca,'fontsize',20);
    
    yyaxis left
    yy_jet_left = pixel_size*(max(ymax_LeftJet(Left_jet_index )) - ymax_LeftJet(Left_jet_index ) );
    yy_jet_right = pixel_size*(max(ymax_RightJet(Right_jet_index)) - ymax_RightJet(Right_jet_index) );
%     yy_jet_mid = (round(rect(4))-ymax_MidJet)*pixel_size;

    plot(T_all_nogap(Left_jet_index )/60, pixel_size*(max(ymax_LeftJet(Left_jet_index )) - ymax_LeftJet(Left_jet_index ) ), 'ro','linewidth',2);
    plot(T_all_nogap(Right_jet_index )/60, pixel_size*(max(ymax_RightJet(Right_jet_index)) - ymax_RightJet(Right_jet_index) ), 'bo','linewidth',2);
%     plot(T_all_nogap(Mid_jet_index )/60, yy_jet_mid(Mid_jet_index )-min(yy_jet_mid(Mid_jet_index) ), 'mo','linewidth',2);

    set(gca,'YColor',[0 0 0]); 
    box on; 

    
%     set(gca, 'YScale', 'log')
%     set(gca, 'XScale', 'log')
    ylabel('y_{top} [mm]')
    xlabel('Time [Hours]')
    
    yyaxis right 
    plot(T_all_nogap(1:gap:nb_all_img)/60, I_ave(1:gap:nb_all_img)/255, 'k+-','linewidth',2);
    ylabel('Average intenstiy')
    set(gca,'YColor',[0 0 0]); 
    
    legend('Left front','Right front','Bacteria intensity','location','northwest');
    legend boxoff 
    
    set(gcf,'PaperPositionMode','Manual')
    set(gcf,'PaperUnits','inches')
    
    set(gcf,'PaperSize',[12 6])
    set(gcf,'PaperPosition',[0 0 12 6])
    saveas(gcf, ['jets_z_total_intensity_normal'], 'pdf')
 
    
 %% plot the height of the cloud as a function of time 
    
    figure; hold on
    set(gca,'fontsize',20);
    
    yyaxis left
    yy_jet_left = pixel_size*(max(ymax_LeftJet(Left_jet_index )) - ymax_LeftJet(Left_jet_index ) );
    yy_jet_right = pixel_size*(max(ymax_RightJet(Right_jet_index)) - ymax_RightJet(Right_jet_index) );
%     yy_jet_mid = (round(rect(4))-ymax_MidJet)*pixel_size;

    plot(T_all_nogap(Left_jet_index )/60, pixel_size*(max(ymax_LeftJet(Left_jet_index )) - ymax_LeftJet(Left_jet_index ) ), 'ro','linewidth',2);
    plot(T_all_nogap(Right_jet_index )/60, pixel_size*(max(ymax_RightJet(Right_jet_index)) - ymax_RightJet(Right_jet_index) ), 'bo','linewidth',2);
%     plot(T_all_nogap(Mid_jet_index )/60, yy_jet_mid(Mid_jet_index )-min(yy_jet_mid(Mid_jet_index) ), 'mo','linewidth',2);

    set(gca,'YColor',[0 0 0]); 
    box on; 

    
    set(gca, 'YScale', 'log')
    set(gca, 'XScale', 'log')
    ylabel('y_{top} [mm]')
    xlabel('Time [Hours]')
    
    yyaxis right 
    plot(T_all_nogap(1:gap:nb_all_img)/60, I_ave(1:gap:nb_all_img)/255, 'k+-','linewidth',2);
    ylabel('Average intenstiy')
    set(gca,'YColor',[0 0 0]); 
    
    legend('Left front','Right front','Bacteria intensity','location','northwest');
    legend boxoff 
    
    set(gcf,'PaperPositionMode','Manual')
    set(gcf,'PaperUnits','inches')
    
    set(gcf,'PaperSize',[12 6])
    set(gcf,'PaperPosition',[0 0 12 6])
    saveas(gcf, ['jets_z_total_intensity'], 'pdf')
 
    
    %% plot the width of the jets as a function of time 
    
    figure; hold on
    set(gca,'fontsize',20);
    

    plot(T_all_nogap(Left_jet_index )/60, width_Ljet_mean(Left_jet_index)*pixel_size, 'ro','linewidth',2);
    plot(T_all_nogap(Right_jet_index )/60, width_Rjet_mean(Right_jet_index)*pixel_size, 'bo','linewidth',2);
%     plot(T_all_nogap(Mid_jet_index )/60, width_Mjet_mean(Mid_jet_index)*pixel_size, 'mo','linewidth',2);

    
    box on; 
    legend('Left jet','Right jet','Middle jet','location','northwest');
    legend boxoff 
    ylabel('Width [mm]')
    xlabel('Time [Hours]')
    set(gcf,'PaperPositionMode','Manual')
    set(gcf,'PaperUnits','inches')
    
    set(gcf,'PaperSize',[12 6])
    set(gcf,'PaperPosition',[0 0 12 6])
    saveas(gcf, ['Jets_width_log_oct'], 'pdf')
    
    
    save('right_chamber2.mat','-v7.3','-nocompression')
