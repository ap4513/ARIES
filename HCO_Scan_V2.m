close all;
clear all;

[file,path] = uigetfile('MultiSelect','on');

% load files
for i = size(file)
    load(fullfile(path,string(file(i))));
end



%% PARAMETERS

% Range for lineout
range = 20;


% Setup window for filtering
windowSize = 10;
d = 1;
c = (1/windowSize)*ones(1,windowSize);

% Set prominence for finding HCO position
prom = 150;

% Moving average parameter for HCO data (Select based on data
movav = 3;

%% DATA ANALYSIS

HCO_locs = zeros(1, size(I_av,3));

    
for i = 1:size(I_av,3)

    image = squeeze(I_av(:,:,i));
    lineout(i,:) = sum(image(140:140+range,:))/(range+1);

    % Filter out the noise and remove background
    lineout_filtered(i,:) =filter(c,d,lineout(i,:));
    [peaks,locations] = findpeaks(lineout_filtered(i,5:end), 'MinPeakProminence',prom);

    if isempty(locations)
      break; 
    end

    HCO_locs(i) = max(locations);

    % Plots working (slows down program)
%         subplot(2,1,1);
%         
%        imagesc(image(:,5:end));
%        axis([0 1000 0 300]);
%        drawnow
%        title(num2str(i));
%              
%        subplot(2,1,2);
%        plot(lineout_filtered(5:end));
%        findpeaks(lineout_filtered(i,5:end), 'MinPeakProminence',100);
%        axis([0 1000 00 3000]);
%        drawnow


end

HCO_Pos = HCO_locs(HCO_locs ~= 0);


%% PLOTTING


% Calculates thickness of wedges
% Thickness = 2.*((0.8 .* (WedgePositions(1:length(HCO_locs))-min(WedgePositions))/200).* tand(2.3)+0.1);
Thickness = 2.*((0.8 .* (WedgePositions(1:length(HCO_locs))-min(WedgePositions))/200).* tand(2.3));

% Calculates CEP from thickness
CEP = Thickness./ (0.0491./2);



figure;
h =imagesc(lineout_filtered');
set(gca,'YDir','normal');
set(h, 'XData', CEP);
axis tight
hold on;
plot(CEP,HCO_locs,'k','linewidth',1);
colorbar;
xlabel('CEP (\pi rads)');
ylabel('Energy (pixels)');
% str = strcat('HCO Positions');
title('HCO Positions');


figure;
h =imagesc(lineout_filtered');
set(gca,'YDir','normal');
set(h, 'XData', CEP);
axis tight
hold on;
plot(CEP,HCO_locs,'k','linewidth',1);
ylim([size(lineout_filtered,2)/2, size(lineout_filtered,2)]);
colorbar;
xlabel('CEP (\pi rads)');
ylabel('Energy (pixels)');
title('HCO Positions (Zoomed)');



HCO_locs_movaverag = movmean(HCO_locs,movav);
figure;
h =imagesc(lineout_filtered');
set(gca,'YDir','normal');
set(h, 'XData', CEP);
axis tight
hold on;
plot(CEP,HCO_locs_movaverag,'k','linewidth',1);
ylim([size(lineout_filtered,2)/2, size(lineout_filtered,2)]);
colorbar;
xlabel('CEP (\pi rads)');
ylabel('Energy (pixels)');
title('HCO Positions (Averaged, Zoomed)');

% Filter Test
% windowSize = 4;
% d = 1;
% c = (1/windowSize)*ones(1,windowSize);
% 
% figure;
% filteredstuff = filter(c,d,HCO_locs);
% h =imagesc(lineout_filtered');
% set(gca,'YDir','normal');
% set(h, 'XData', CEP);
% axis tight
% hold on;
% plot(CEP,filteredstuff,'k','linewidth',1);
% ylim([size(lineout_filtered,2)/2, size(lineout_filtered,2)]);
% colorbar;
% xlabel('CEP (\pi rads)');
% ylabel('Energy (pixels)');
% 

























