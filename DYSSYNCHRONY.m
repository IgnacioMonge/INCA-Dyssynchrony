%% LV MECHANICAL DYSSYNCHRONY ANALYZER FOR INCA SYSTEM
%% BY MANUEL IGNACIO MONGE GARCÍA - ignacio_monge@edwards.com

clearvars; close all; clc;

Fs = 250;
dt = 1/Fs;
negative_values = 0;

%% loading interface for LV data from the INCA System
% Please, adjust according to the file exported from INCA software

[filename,path] = uigetfile('*.csv');
data_INCA= dlmread([path,filename], ',', 3, 0);

    beat_n = data_INCA (:,1)';
    beat_time = data_INCA (:,2)';
    Segment1_INCA = data_INCA (:,3);
    Segment2_INCA = data_INCA (:,4);
    Segment3_INCA = data_INCA (:,5);
    Segment4_INCA = data_INCA (:,6);
    Segment5_INCA = data_INCA (:,7);
    ECG_INCA = data_INCA (:,8)';
    pressure_INCA = data_INCA (:,9)';
    volume_INCA = data_INCA (:,15)';
    dPdt_INCA = data_INCA(:,16);

% Constructing vector of interpolation points
DT1 = length(pressure_INCA);
tt1_interp = linspace(0,DT1*dt,DT1);

% Determining the threshold for R peak finder
 [pks_beats,loc_beats]=findpeaks(volume_INCA,'MinPeakProminence',std(volume_INCA)); 

 % Finding mean cycles duration
delta1 = zeros(1,length(loc_beats)-1);
for kk1 = 1:(length(loc_beats)-1)
    delta1(kk1) = (loc_beats(kk1+1)) - loc_beats(kk1);
end
 time_threshold = round(min(rmoutliers(delta1,'median')));
 
%% EKG Analsis
ecg_max = max(rmoutliers(ECG_INCA, 'gesd'));
thresholdx = ecg_max*0.5;
[rpeaks,locs_Rpeaks] = findpeaks(ECG_INCA,'MaxPeakWidth',50,'MinPeakHeight',thresholdx,'MinPeakDistance',time_threshold*0.8);
Rwaves =tt1_interp(locs_Rpeaks);
figure (1001);plot(ECG_INCA,'DisplayName','ecg_tonom');hold on;plot(locs_Rpeaks,rpeaks,'*');hold off

% Finding mean cycles duration
delta2 = zeros(1,length(locs_Rpeaks)-1);
for kk1 = 1:(length(locs_Rpeaks)-1)
    delta2(kk1) = (locs_Rpeaks(kk1+1)) - locs_Rpeaks(kk1);
end
min_duration = round(min(rmoutliers(delta2,'gesd')));

%% Time varying elastance calculation
Vd = min(volume_INCA) - (0.5*(max(volume_INCA)-min(volume_INCA)));
LV_elastance = pressure_INCA./(volume_INCA-Vd);

%% Max LV Pressure
max_prominence_P = median(pressure_INCA);
[pks_pressure,locs_pressure]=findpeaks(pressure_INCA,'MinPeakDistance',min_duration*0.8);
systolic_pressure = pks_pressure;
t_systolicpressure = tt1_interp(locs_pressure);

%% Determine EDP and EDV at R wave on the EKG
t_end_diastole = Rwaves;
ED_volume = interp1(tt1_interp,volume_INCA,Rwaves);
ED_pressure = interp1(tt1_interp,pressure_INCA,Rwaves);

% Determine the Emax and t_Emax
 max_prominence = median(LV_elastance);
 [pks,locs]=findpeaks(LV_elastance,'MinPeakDistance',min_duration*0.8);
 emax = pks;
 time_emax = tt1_interp(locs);
 pressure_emax = pressure_INCA(locs);
 volume_emax = volume_INCA(locs);
 
 % Determine LV pressure and volume at t_Emax
 LV_ESP = interp1(tt1_interp,pressure_INCA,time_emax);
 LV_ESV = interp1(tt1_interp,volume_INCA,time_emax);
 
%  Correct the offsset of ED and ES values 
 if length(ED_volume) > length(LV_ESP)
     if t_end_diastole(end) > time_emax(end)
         t_end_diastole = t_end_diastole(1:end-1);
         ED_volume = ED_volume(1:end-1);
         ED_pressure = ED_pressure(1:end-1);
     end
 end
  
%& Derivative analyss for dyssynchrony analysis
 LV_dVdt_INCA = diff(volume_INCA);
 Seg1_dVdt_INCA = gradient(Segment1_INCA);
 Seg2_dVdt_INCA = gradient(Segment2_INCA);
 Seg3_dVdt_INCA = gradient(Segment3_INCA);
 Seg4_dVdt_INCA = gradient(Segment4_INCA);
 Seg5_dVdt_INCA = gradient(Segment5_INCA);
%  Seg6_dVdt_INCA = gradient(Segment6_INCA,tt1_interp);

% Total LV volume
LV_vol_segments = [Segment1_INCA Segment2_INCA Segment3_INCA Segment4_INCA Segment5_INCA]';
total_LVvolume = sum(LV_vol_segments);
dVdt_LVtotal = gradient(total_LVvolume);

% Dyssynchrony Total LV Volume vs Segmental Volume
Seg1_discrepancy = Seg1_dVdt_INCA' .*dVdt_LVtotal;
Seg2_discrepancy = Seg2_dVdt_INCA' .*dVdt_LVtotal;
Seg3_discrepancy = Seg3_dVdt_INCA' .*dVdt_LVtotal;
Seg4_discrepancy = Seg4_dVdt_INCA' .*dVdt_LVtotal;
Seg5_discrepancy = Seg5_dVdt_INCA' .*dVdt_LVtotal;

%& Segmental IF for Internal Flow analysis
Segment_dVdt = [abs(Seg1_dVdt_INCA) abs(Seg2_dVdt_INCA) abs(Seg3_dVdt_INCA) abs(Seg4_dVdt_INCA) abs(Seg5_dVdt_INCA)];
Segment_dVdt  = Segment_dVdt';
A = sum(Segment_dVdt);
B = abs(dVdt_LVtotal);
INTERNAL_FLOW = (A-B)/2;

% Segmental Internal Flow
Segm1_IF = (abs(Seg1_dVdt_INCA)'-B)/2;
Segm2_IF = (abs(Seg2_dVdt_INCA)'-B)/2;
Segm3_IF = (abs(Seg3_dVdt_INCA)'-B)/2;
Segm4_IF = (abs(Seg4_dVdt_INCA)'-B)/2;
Segm5_IF = (abs(Seg5_dVdt_INCA)'-B)/2;


%% CYCLES SEGMENTATION 
n_cycles = length(t_end_diastole)-1;
 
volumes = zeros(1,1000);
for cycle = 1:length(t_end_diastole)-1
    volume_per_cycle = volume_INCA(tt1_interp<t_end_diastole(cycle+1) & tt1_interp>t_end_diastole(cycle));
    volumes (cycle,1:length(volume_per_cycle)) = volume_per_cycle;
end

% Segmental volume analysis
segmt1 = zeros(1,2000);
for cycle = 1:length(t_end_diastole)-1
    segmt1_per_cycle = Segment1_INCA(tt1_interp<t_end_diastole(cycle+1) & tt1_interp>t_end_diastole(cycle));
    segmt1 (cycle,1:length(segmt1_per_cycle)) = segmt1_per_cycle;
end
segmt2 = zeros(1,2000);
for cycle = 1:length(t_end_diastole)-1
    segmt2_per_cycle = Segment2_INCA(tt1_interp<t_end_diastole(cycle+1) & tt1_interp>t_end_diastole(cycle));
    segmt2 (cycle,1:length(segmt2_per_cycle)) = segmt2_per_cycle;
end
segmt3 = zeros(1,2000);
for cycle = 1:length(t_end_diastole)-1
    segmt3_per_cycle = Segment3_INCA(tt1_interp<t_end_diastole(cycle+1) & tt1_interp>t_end_diastole(cycle));
    segmt3 (cycle,1:length(segmt3_per_cycle)) = segmt3_per_cycle;
end
segmt4 = zeros(1,2000);
for cycle = 1:length(t_end_diastole)-1
    segmt4_per_cycle = Segment4_INCA(tt1_interp<t_end_diastole(cycle+1) & tt1_interp>t_end_diastole(cycle));
    segmt4 (cycle,1:length(segmt4_per_cycle)) = segmt4_per_cycle;
end
segmt5 = zeros(1,2000);
for cycle = 1:length(t_end_diastole)-1
    segmt5_per_cycle = Segment5_INCA(tt1_interp<t_end_diastole(cycle+1) & tt1_interp>t_end_diastole(cycle));
    segmt5 (cycle,1:length(segmt5_per_cycle)) = segmt5_per_cycle;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Dyssynchrony analysis

DYS_segmt1 = zeros(1,2000);
for cycle = 1:length(t_end_diastole)-1
    DYS_segmt1_per_cycle = Seg1_discrepancy(tt1_interp<t_end_diastole(cycle+1) & tt1_interp>t_end_diastole(cycle));
    DYS_segmt1 (cycle,1:length(DYS_segmt1_per_cycle )) = DYS_segmt1_per_cycle ;
end
DYS_segmt2 = zeros(1,2000);
for cycle = 1:length(t_end_diastole)-1
    DYS_segmt2_per_cycle = Seg2_discrepancy(tt1_interp<t_end_diastole(cycle+1) & tt1_interp>t_end_diastole(cycle));
    DYS_segmt2 (cycle,1:length(DYS_segmt2_per_cycle)) = DYS_segmt2_per_cycle;
end
DYS_segmt3 = zeros(1,2000);
for cycle = 1:length(t_end_diastole)-1
    DYS_segmt3_per_cycle = Seg3_discrepancy(tt1_interp<t_end_diastole(cycle+1) & tt1_interp>t_end_diastole(cycle));
    DYS_segmt3 (cycle,1:length(DYS_segmt3_per_cycle)) = DYS_segmt3_per_cycle;
end
DYS_segmt4 = zeros(1,2000);
for cycle = 1:length(t_end_diastole)-1
    DYS_segmt4_per_cycle = Seg4_discrepancy(tt1_interp<t_end_diastole(cycle+1) & tt1_interp>t_end_diastole(cycle));
    DYS_segmt4 (cycle,1:length(DYS_segmt4_per_cycle)) = DYS_segmt4_per_cycle;
end
DYS_segmt5 = zeros(1,2000);
for cycle = 1:length(t_end_diastole)-1
    DYS_segmt5_per_cycle = Seg5_discrepancy(tt1_interp<t_end_diastole(cycle+1) & tt1_interp>t_end_diastole(cycle));
    DYS_segmt5 (cycle,1:length(DYS_segmt5_per_cycle)) = DYS_segmt5_per_cycle;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INTERNAL FLOW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IFs = zeros(1,1000);
for cycle = 1:n_cycles
    IF_per_cycle = INTERNAL_FLOW(tt1_interp<t_end_diastole(cycle+1) & tt1_interp>t_end_diastole(cycle));
    IFs(cycle,1:length(IF_per_cycle)) = IF_per_cycle;
end

IFF_segmt1 = zeros(1,2000);
for cycle = 1:n_cycles
    IFF_segmt1_per_cycle = Segm1_IF(tt1_interp<t_end_diastole(cycle+1) & tt1_interp>t_end_diastole(cycle));
    IFF_segmt1 (cycle,1:length(IFF_segmt1_per_cycle )) = IFF_segmt1_per_cycle ;
end
IFF_segmt2 = zeros(1,2000);
for cycle = 1:n_cycles
    IFF_segmt2_per_cycle = Segm2_IF(tt1_interp<t_end_diastole(cycle+1) & tt1_interp>t_end_diastole(cycle));
    IFF_segmt2 (cycle,1:length(IFF_segmt2_per_cycle )) = IFF_segmt2_per_cycle ;
end
IFF_segmt3 = zeros(1,2000);
for cycle = 1:n_cycles
    IFF_segmt3_per_cycle = Segm3_IF(tt1_interp<t_end_diastole(cycle+1) & tt1_interp>t_end_diastole(cycle));
    IFF_segmt3 (cycle,1:length(IFF_segmt3_per_cycle )) = IFF_segmt3_per_cycle ;
end
IFF_segmt4 = zeros(1,2000);
for cycle = 1:n_cycles
    IFF_segmt4_per_cycle = Segm4_IF(tt1_interp<t_end_diastole(cycle+1) & tt1_interp>t_end_diastole(cycle));
    IFF_segmt4 (cycle,1:length(IFF_segmt4_per_cycle )) = IFF_segmt4_per_cycle ;
end
IFF_segmt5 = zeros(1,2000);
for cycle = 1:n_cycles
    IFF_segmt5_per_cycle = Segm5_IF(tt1_interp<t_end_diastole(cycle+1) & tt1_interp>t_end_diastole(cycle));
    IFF_segmt5 (cycle,1:length(IFF_segmt5_per_cycle)) = IFF_segmt5_per_cycle;
end

pressures = zeros(1,1000);
for cycle = 1:length(t_end_diastole)-1
    pressure_per_cycle = pressure_INCA(tt1_interp<t_end_diastole(cycle+1) & tt1_interp>t_end_diastole(cycle));
    pressures (cycle,1:length(pressure_per_cycle)) = pressure_per_cycle;
end

elastances = zeros(1,1000);
for cycle = 1:length(t_end_diastole)-1
    elastance_per_cycle = LV_elastance(tt1_interp<t_end_diastole(cycle+1) & tt1_interp>t_end_diastole(cycle));
    elastances (cycle,1:length(elastance_per_cycle)) = elastance_per_cycle;
end

ecg = zeros(1,1000);
for cycle = 1:length(t_end_diastole)-1
    ecg_per_cycle = ECG_INCA (tt1_interp<t_end_diastole(cycle+1) & tt1_interp>t_end_diastole(cycle));
    ecg (cycle,1:length(ecg_per_cycle)) = ecg_per_cycle;
end

dPdts = zeros(1,1000);
for cycle = 1:length(t_end_diastole)-1
    dPdt_per_cycle = dPdt_INCA(tt1_interp<t_end_diastole(cycle+1) & tt1_interp>t_end_diastole(cycle));
    dPdts(cycle,1:length(dPdt_per_cycle)) = dPdt_per_cycle;
end

dVdts = zeros(1,1000);
for cycle = 1:length(t_end_diastole)-1
    dVdt_per_cycle = LV_dVdt_INCA(tt1_interp<t_end_diastole(cycle+1) & tt1_interp>t_end_diastole(cycle));
    dVdts(cycle,1:length(dVdt_per_cycle)) = dVdt_per_cycle;
end

dVdt_uncals = zeros(1,1000);
for cycle = 1:length(t_end_diastole)-1
    dVdt_uncal_per_cycle = dVdt_LVtotal(tt1_interp<t_end_diastole(cycle+1) & tt1_interp>t_end_diastole(cycle));
    dVdt_uncals(cycle,1:length(dVdt_uncal_per_cycle)) = dVdt_uncal_per_cycle;
end

timelines = zeros(1,1000);
for cycle = 1:length(t_end_diastole)-1
    time_per_cycle = tt1_interp(tt1_interp<t_end_diastole(cycle+1) & tt1_interp>t_end_diastole(cycle));
    timelines (cycle,1:length(time_per_cycle)) = time_per_cycle;
end

SW_cycles = zeros(1,length(cycle));
EDV_cycles  = zeros(1,length(cycle)); 
EDP_cycles  = zeros(1,length(cycle));
ESV_cycles  = zeros(1,length(cycle)); 
ESP_cycles  = zeros(1,length(cycle));
tED_cycles = zeros(1,length(cycle)); 
max_elastance = zeros(1,length(cycle));
t_Ees_max = zeros(1,length(cycle));

IFF_cycles = zeros(1,length(cycle));
IFFsystole_cycle = zeros(1,length(cycle));
IFFdiastole_cycle = zeros(1,length(cycle));
IFF_Segment1_cycles = zeros(1,length(cycle));
IFF_Segment2_cycles = zeros(1,length(cycle));
IFF_Segment3_cycles = zeros(1,length(cycle));
IFF_Segment4_cycles = zeros(1,length(cycle));
IFF_Segment5_cycles = zeros(1,length(cycle));

DYS_Seg1_cycles = zeros(1,length(cycle));
DYS_Seg2_cycles = zeros(1,length(cycle));
DYS_Seg3_cycles = zeros(1,length(cycle));
DYS_Seg4_cycles = zeros(1,length(cycle));
DYS_Seg5_cycles = zeros(1,length(cycle));

DYStotal_seg1 = zeros(1,length(cycle));
DYStotal_seg2 = zeros(1,length(cycle));
DYStotal_seg3 = zeros(1,length(cycle));
DYStotal_seg4 = zeros(1,length(cycle));
DYStotal_seg5 = zeros(1,length(cycle));

DYSsystolic_seg1 = zeros(1,length(cycle));
DYSsystolic_seg2 = zeros(1,length(cycle));
DYSsystolic_seg3 = zeros(1,length(cycle));
DYSsystolic_seg4 = zeros(1,length(cycle));
DYSsystolic_seg5 = zeros(1,length(cycle));

DYSdiastolic_seg1 = zeros(1,length(cycle));
DYSdiastolic_seg2 = zeros(1,length(cycle));
DYSdiastolic_seg3 = zeros(1,length(cycle));
DYSdiastolic_seg4 = zeros(1,length(cycle));
DYSdiastolic_seg5 = zeros(1,length(cycle));

DYSearly_systolic_seg1 = zeros(1,length(cycle));
DYSearly_systolic_seg2 = zeros(1,length(cycle));
DYSearly_systolic_seg3 = zeros(1,length(cycle));
DYSearly_systolic_seg4 = zeros(1,length(cycle));
DYSearly_systolic_seg5 = zeros(1,length(cycle));

DYSlate_systolic_seg1 = zeros(1,length(cycle));
DYSlate_systolic_seg2 = zeros(1,length(cycle));
DYSlate_systolic_seg3 = zeros(1,length(cycle));
DYSlate_systolic_seg4 = zeros(1,length(cycle));
DYSlate_systolic_seg5 = zeros(1,length(cycle));

DYSearly_diastolic_seg1 = zeros(1,length(cycle));
DYSearly_diastolic_seg2 = zeros(1,length(cycle));
DYSearly_diastolic_seg3 = zeros(1,length(cycle));
DYSearly_diastolic_seg4 = zeros(1,length(cycle));
DYSearly_diastolic_seg5 = zeros(1,length(cycle));

DYSlate_diastolic_seg1 = zeros(1,length(cycle));
DYSlate_diastolic_seg2 = zeros(1,length(cycle));
DYSlate_diastolic_seg3 = zeros(1,length(cycle));
DYSlate_diastolic_seg4 = zeros(1,length(cycle));
DYSlate_diastolic_seg5 = zeros(1,length(cycle));

DYS_TOTAL = zeros(1,length(cycle));
DYS_SYSTOLIC = zeros(1,length(cycle));
DYS_DIASTOLIC = zeros(1,length(cycle));
DYS_SYSearly = zeros(1,length(cycle));
DYS_SYSlate = zeros(1,length(cycle));
DYS_DIAearly = zeros(1,length(cycle));
DYS_DIAlate = zeros(1,length(cycle));

ENTROPY_ECG = zeros(1,length(cycle));
SDI_cycle = zeros (1,length(cycle));
SDInew_cycle = zeros (1,length(cycle));

%% BEAT TO BEAT ANALYSIS

for cycle = 1:length(t_end_diastole)-1
volume_cycles = volumes(cycle,:);
volume_cycles = volume_cycles(volume_cycles~=0);

%segmental analysis
segm1_cycles = segmt1(cycle,:);
segm1_cycles = segm1_cycles(segm1_cycles~=0);
segm2_cycles = segmt2(cycle,:);
segm2_cycles = segm2_cycles(segm2_cycles~=0);
segm3_cycles = segmt3(cycle,:);
segm3_cycles = segm3_cycles(segm3_cycles~=0);
segm4_cycles = segmt4(cycle,:);
segm4_cycles = segm4_cycles(segm4_cycles~=0);
segm5_cycles = segmt5(cycle,:);
segm5_cycles = segm5_cycles(segm5_cycles~=0);

%dyssynchrony analysis
DYS_Seg1_cycles = DYS_segmt1(cycle,:);
DYS_Seg1_cycles = DYS_Seg1_cycles(DYS_Seg1_cycles~=0);
DYS_Seg2_cycles = DYS_segmt2(cycle,:);
DYS_Seg2_cycles = DYS_Seg2_cycles(DYS_Seg2_cycles~=0);
DYS_Seg3_cycles = DYS_segmt3(cycle,:);
DYS_Seg3_cycles = DYS_Seg3_cycles(DYS_Seg3_cycles~=0);
DYS_Seg4_cycles = DYS_segmt4(cycle,:);
DYS_Seg4_cycles = DYS_Seg4_cycles(DYS_Seg4_cycles~=0);
DYS_Seg5_cycles = DYS_segmt5(cycle,:);
DYS_Seg5_cycles = DYS_Seg5_cycles(DYS_Seg5_cycles~=0);

dVdt_per_cycle = dVdts(cycle,:);
dVdt_per_cycle = dVdt_per_cycle(dVdt_per_cycle~=0);
pressure_cycles = pressures(cycle,:); % -- use volume instead pressure 
pressure_cycles = pressure_cycles(volume_cycles~=0); % we do exclude any actual pressure = 0 value
dPdt_cycles = dPdts(cycle,:);
dPdt_cycles = dPdt_cycles(dPdt_cycles~=0);
elastances_cycles = elastances(cycle,:);
elastances_cycles = elastances_cycles(volume_cycles~=0);
timeline_cycle = timelines(cycle,:);
timeline_cycles = timeline_cycle(timeline_cycle~=0)';
ecg_cycle = ecg(cycle,:);
ecg_cycles = ecg_cycle(ecg_cycle~=0)';
IF_cycle = IFs(cycle,:);
IF_cycles = IF_cycle(IF_cycle~=0)';
dVdt_uncal_cycles = dVdt_uncals(cycle,:);
dVdt_uncal_cycles = dVdt_uncal_cycles(dVdt_uncal_cycles~=0)';

EDV_cycles (cycle) = volume_cycles(end);
EDP_cycles (cycle) = pressure_cycles(end);
tED_cycles (cycle) = timeline_cycles(end);
[max_elastance(cycle),Emax_cycles_idx] = max(elastances_cycles);
t_Ees_max (cycle) = timeline_cycles(Emax_cycles_idx);
[~,dPdtmax_idx] = max(dPdt_cycles);
[~,dPdtmin_idx] = min(dPdt_cycles);
ESP_cycles (cycle) = pressure_cycles(Emax_cycles_idx);
ESV_cycles (cycle) = volume_cycles(Emax_cycles_idx);

%%DYSYNCHRONY ANALYSIS
% Total Segmental Dyssynchrony
DYStotal_seg1 (cycle) = (length(find(DYS_Seg1_cycles<0))/length(DYS_Seg1_cycles))*100;
DYStotal_seg2 (cycle) = (length(find(DYS_Seg2_cycles<0))/length(DYS_Seg2_cycles))*100;
DYStotal_seg3 (cycle) = (length(find(DYS_Seg3_cycles<0))/length(DYS_Seg3_cycles))*100;
DYStotal_seg4 (cycle) = (length(find(DYS_Seg4_cycles<0))/length(DYS_Seg4_cycles))*100;
DYStotal_seg5 (cycle) = (length(find(DYS_Seg5_cycles<0))/length(DYS_Seg5_cycles))*100;

% Systolic Segmental Dyssynchrony
Systolic_DYS_Seg1_cycles = DYS_Seg1_cycles(1:Emax_cycles_idx);
Systolic_DYS_Seg2_cycles = DYS_Seg2_cycles(1:Emax_cycles_idx);
Systolic_DYS_Seg3_cycles = DYS_Seg3_cycles(1:Emax_cycles_idx);
Systolic_DYS_Seg4_cycles = DYS_Seg4_cycles(1:Emax_cycles_idx);
Systolic_DYS_Seg5_cycles = DYS_Seg5_cycles(1:Emax_cycles_idx);

DYSsystolic_seg1 (cycle) = (length(find(Systolic_DYS_Seg1_cycles<0))/length(Systolic_DYS_Seg1_cycles))*100;
DYSsystolic_seg2 (cycle) = (length(find(Systolic_DYS_Seg2_cycles<0))/length(Systolic_DYS_Seg2_cycles))*100;
DYSsystolic_seg3 (cycle) = (length(find(Systolic_DYS_Seg3_cycles<0))/length(Systolic_DYS_Seg3_cycles))*100;
DYSsystolic_seg4 (cycle) = (length(find(Systolic_DYS_Seg4_cycles<0))/length(Systolic_DYS_Seg4_cycles))*100;
DYSsystolic_seg5 (cycle) = (length(find(Systolic_DYS_Seg5_cycles<0))/length(Systolic_DYS_Seg5_cycles))*100;

% Diastolic Segmental Dyssynchrony
Diastolic_DYS_Seg1_cycles = DYS_Seg1_cycles(Emax_cycles_idx:end);
Diastolic_DYS_Seg2_cycles = DYS_Seg2_cycles(Emax_cycles_idx:end);
Diastolic_DYS_Seg3_cycles = DYS_Seg3_cycles(Emax_cycles_idx:end);
Diastolic_DYS_Seg4_cycles = DYS_Seg4_cycles(Emax_cycles_idx:end);
Diastolic_DYS_Seg5_cycles = DYS_Seg5_cycles(Emax_cycles_idx:end);

DYSdiastolic_seg1 (cycle) = (length(find(Diastolic_DYS_Seg1_cycles<0))/length(Diastolic_DYS_Seg1_cycles))*100;
DYSdiastolic_seg2 (cycle) = (length(find(Diastolic_DYS_Seg2_cycles<0))/length(Diastolic_DYS_Seg2_cycles))*100;
DYSdiastolic_seg3 (cycle) = (length(find(Diastolic_DYS_Seg3_cycles<0))/length(Diastolic_DYS_Seg3_cycles))*100;
DYSdiastolic_seg4 (cycle) = (length(find(Diastolic_DYS_Seg4_cycles<0))/length(Diastolic_DYS_Seg4_cycles))*100;
DYSdiastolic_seg5 (cycle) = (length(find(Diastolic_DYS_Seg5_cycles<0))/length(Diastolic_DYS_Seg5_cycles))*100;

% Early Segmental Dyssynchrony (from EDP to dPdtmax)
Early_Systolic_DYS_Seg1_cycles = DYS_Seg1_cycles(1:dPdtmax_idx);
Early_Systolic_DYS_Seg2_cycles = DYS_Seg2_cycles(1:dPdtmax_idx);
Early_Systolic_DYS_Seg3_cycles = DYS_Seg3_cycles(1:dPdtmax_idx);
Early_Systolic_DYS_Seg4_cycles = DYS_Seg4_cycles(1:dPdtmax_idx);
Early_Systolic_DYS_Seg5_cycles = DYS_Seg5_cycles(1:dPdtmax_idx);

DYSearly_systolic_seg1 (cycle) = (length(find(Early_Systolic_DYS_Seg1_cycles<0))/length(Early_Systolic_DYS_Seg1_cycles))*100;
DYSearly_systolic_seg2 (cycle) = (length(find(Early_Systolic_DYS_Seg2_cycles<0))/length(Early_Systolic_DYS_Seg2_cycles))*100;
DYSearly_systolic_seg3 (cycle) = (length(find(Early_Systolic_DYS_Seg3_cycles<0))/length(Early_Systolic_DYS_Seg3_cycles))*100;
DYSearly_systolic_seg4 (cycle) = (length(find(Early_Systolic_DYS_Seg4_cycles<0))/length(Early_Systolic_DYS_Seg4_cycles))*100;
DYSearly_systolic_seg5 (cycle) = (length(find(Early_Systolic_DYS_Seg5_cycles<0))/length(Early_Systolic_DYS_Seg5_cycles))*100;

% Early Segmental Dyssynchrony (from EDP to dPdtmax)
Late_Systolic_DYS_Seg1_cycles = DYS_Seg1_cycles(dPdtmax_idx:Emax_cycles_idx);
Late_Systolic_DYS_Seg2_cycles = DYS_Seg2_cycles(dPdtmax_idx:Emax_cycles_idx);
Late_Systolic_DYS_Seg3_cycles = DYS_Seg3_cycles(dPdtmax_idx:Emax_cycles_idx);
Late_Systolic_DYS_Seg4_cycles = DYS_Seg4_cycles(dPdtmax_idx:Emax_cycles_idx);
Late_Systolic_DYS_Seg5_cycles = DYS_Seg5_cycles(dPdtmax_idx:Emax_cycles_idx);

DYSlate_systolic_seg1 (cycle) = (length(find(Late_Systolic_DYS_Seg1_cycles<0))/length(Late_Systolic_DYS_Seg1_cycles))*100;
DYSlate_systolic_seg2 (cycle) = (length(find(Late_Systolic_DYS_Seg2_cycles<0))/length(Late_Systolic_DYS_Seg2_cycles))*100;
DYSlate_systolic_seg3 (cycle) = (length(find(Late_Systolic_DYS_Seg3_cycles<0))/length(Late_Systolic_DYS_Seg3_cycles))*100;
DYSlate_systolic_seg4 (cycle) = (length(find(Late_Systolic_DYS_Seg4_cycles<0))/length(Late_Systolic_DYS_Seg4_cycles))*100;
DYSlate_systolic_seg5 (cycle) = (length(find(Late_Systolic_DYS_Seg5_cycles<0))/length(Late_Systolic_DYS_Seg5_cycles))*100;

% Early Systole Segmental Dyssynchrony DIASTOLE (from EMAX to dPdtmin)
Early_Diastolic_DYS_Seg1_cycles = DYS_Seg1_cycles(Emax_cycles_idx:dPdtmin_idx);
Early_Diastolic_DYS_Seg2_cycles = DYS_Seg2_cycles(Emax_cycles_idx:dPdtmin_idx);
Early_Diastolic_DYS_Seg3_cycles = DYS_Seg3_cycles(Emax_cycles_idx:dPdtmin_idx);
Early_Diastolic_DYS_Seg4_cycles = DYS_Seg4_cycles(Emax_cycles_idx:dPdtmin_idx);
Early_Diastolic_DYS_Seg5_cycles = DYS_Seg5_cycles(Emax_cycles_idx:dPdtmin_idx);

DYSearly_diastolic_seg1 (cycle) = (length(find(Early_Diastolic_DYS_Seg1_cycles<0))/length(Early_Diastolic_DYS_Seg1_cycles))*100;
DYSearly_diastolic_seg2 (cycle) = (length(find(Early_Diastolic_DYS_Seg2_cycles<0))/length(Early_Diastolic_DYS_Seg2_cycles))*100;
DYSearly_diastolic_seg3 (cycle) = (length(find(Early_Diastolic_DYS_Seg3_cycles<0))/length(Early_Diastolic_DYS_Seg3_cycles))*100;
DYSearly_diastolic_seg4 (cycle) = (length(find(Early_Diastolic_DYS_Seg4_cycles<0))/length(Early_Diastolic_DYS_Seg4_cycles))*100;
DYSearly_diastolic_seg5 (cycle) = (length(find(Early_Diastolic_DYS_Seg5_cycles<0))/length(Early_Diastolic_DYS_Seg5_cycles))*100;

% Late Systole Segmental Dyssynchrony DIASTOLE (from dPdtmin to end)
Late_Diastolic_DYS_Seg1_cycles = DYS_Seg1_cycles(dPdtmin_idx:end);
Late_Diastolic_DYS_Seg2_cycles = DYS_Seg2_cycles(dPdtmin_idx:end);
Late_Diastolic_DYS_Seg3_cycles = DYS_Seg3_cycles(dPdtmin_idx:end);
Late_Diastolic_DYS_Seg4_cycles = DYS_Seg4_cycles(dPdtmin_idx:end);
Late_Diastolic_DYS_Seg5_cycles = DYS_Seg5_cycles(dPdtmin_idx:end);

DYSlate_diastolic_seg1 (cycle) = (length(find(Late_Diastolic_DYS_Seg1_cycles<0))/length(Late_Diastolic_DYS_Seg1_cycles))*100;
DYSlate_diastolic_seg2 (cycle) = (length(find(Late_Diastolic_DYS_Seg2_cycles<0))/length(Late_Diastolic_DYS_Seg2_cycles))*100;
DYSlate_diastolic_seg3 (cycle) = (length(find(Late_Diastolic_DYS_Seg3_cycles<0))/length(Late_Diastolic_DYS_Seg3_cycles))*100;
DYSlate_diastolic_seg4 (cycle) = (length(find(Late_Diastolic_DYS_Seg4_cycles<0))/length(Late_Diastolic_DYS_Seg4_cycles))*100;
DYSlate_diastolic_seg5 (cycle) = (length(find(Late_Diastolic_DYS_Seg5_cycles<0))/length(Late_Diastolic_DYS_Seg5_cycles))*100;

% Effective Flow and IF Fraction
Effective_Flow = trapz(abs(dVdt_uncal_cycles));
Systolic_effective_flow = trapz(abs(dVdt_uncal_cycles(1:Emax_cycles_idx)));
Diastolic_effective_flow = trapz(abs(dVdt_uncal_cycles(Emax_cycles_idx:end)));
IFF_cycles (cycle) = (trapz(IF_cycles)/Effective_Flow)*100;
IFFsystole_cycle (cycle) = (trapz(IF_cycles(1:Emax_cycles_idx))/Systolic_effective_flow)*100;
IFFdiastole_cycle (cycle) = (trapz(IF_cycles(Emax_cycles_idx:end))/Diastolic_effective_flow)*100;
IFF_Seg1_cycles = IFF_segmt1(cycle,:);
IFF_Seg1_cycles = IFF_Seg1_cycles(IFF_Seg1_cycles~=0);
IFF_Seg2_cycles = IFF_segmt2(cycle,:);
IFF_Seg2_cycles = IFF_Seg2_cycles(IFF_Seg2_cycles~=0);
IFF_Seg3_cycles = IFF_segmt3(cycle,:);
IFF_Seg3_cycles = IFF_Seg3_cycles(IFF_Seg3_cycles~=0);
IFF_Seg4_cycles = IFF_segmt4(cycle,:);
IFF_Seg4_cycles = IFF_Seg4_cycles(IFF_Seg4_cycles~=0);
IFF_Seg5_cycles = IFF_segmt5(cycle,:);
IFF_Seg5_cycles = IFF_Seg5_cycles(IFF_Seg5_cycles~=0);

% Average calculation of dyssynchrony parameters
% Internal Flow
IFF = mean(IFF_cycles);
IFF_systole = mean(IFFsystole_cycle);
IFF_diastole = mean(IFFdiastole_cycle);
IFF_Segment1 = mean(IFF_Seg1_cycles);
IFF_Segment2 = mean(IFF_Seg2_cycles);
IFF_Segment3 = mean(IFF_Seg3_cycles);
IFF_Segment4 = mean(IFF_Seg4_cycles);
IFF_Segment5 = mean(IFF_Seg5_cycles);

% Total dyssynchrony
DYS_TOTAL (cycle) = (DYStotal_seg1(cycle) +DYStotal_seg2(cycle) + DYStotal_seg3(cycle)  + DYStotal_seg4(cycle) + DYStotal_seg5(cycle)) / 5;
total_DYS = mean(DYS_TOTAL); % Overall dyssychnrony
Segm1_DYStotal = mean(DYStotal_seg1);
Segm2_DYStotal = mean(DYStotal_seg2);
Segm3_DYStotal = mean(DYStotal_seg3);
Segm4_DYStotal = mean(DYStotal_seg4);
Segm5_DYStotal = mean(DYStotal_seg5);

% Systolic dyssynchrony
DYS_SYSTOLIC (cycle) = (DYSsystolic_seg1(cycle) + DYSsystolic_seg2(cycle) + DYSsystolic_seg3(cycle) + DYSsystolic_seg4(cycle) + DYSsystolic_seg5(cycle))/5;
systolic_DYS = mean(DYS_SYSTOLIC); % Overall systolic dyssynchrony
% Segmental analysis
Segm1_DYSsys = mean(DYSsystolic_seg1);
Segm2_DYSsys = mean(DYSsystolic_seg2);
Segm3_DYSsys = mean(DYSsystolic_seg3);
Segm4_DYSsys = mean(DYSsystolic_seg4);
Segm5_DYSsys = mean(DYSsystolic_seg5);
% Early/Late systolic dyssynchrony
DYS_SYSearly (cycle) = (DYSearly_systolic_seg1(cycle) + DYSearly_systolic_seg2(cycle) + DYSearly_systolic_seg3(cycle) + DYSearly_systolic_seg4(cycle) + DYSearly_systolic_seg5(cycle))/5;
early_systole_DYS = mean(DYS_SYSearly);
DYS_SYSlate (cycle) = (DYSlate_systolic_seg1(cycle) + DYSlate_systolic_seg2(cycle) + DYSlate_systolic_seg3(cycle) + DYSlate_systolic_seg4(cycle) + DYSlate_systolic_seg5(cycle))/5;
late_systole_DYS = mean(DYS_SYSlate);

% Diastolic dyssynchrony
DYS_DIASTOLIC (cycle) = (DYSdiastolic_seg1(cycle) + DYSdiastolic_seg2(cycle) + DYSdiastolic_seg3(cycle) + DYSdiastolic_seg4(cycle) + DYSdiastolic_seg5(cycle))/5;
diastolic_DYS = mean(DYS_DIASTOLIC); % Overall diastolic dyssynchrony
% Segmental analysis
Segm1_DYSdia = mean(DYSdiastolic_seg1);
Segm2_DYSdia = mean(DYSdiastolic_seg2);
Segm3_DYSdia = mean(DYSdiastolic_seg3);
Segm4_DYSdia = mean(DYSdiastolic_seg4);
Segm5_DYSdia = mean(DYSdiastolic_seg5);

% Early/Late diastolic dyssynchrony
DYS_DIAearly (cycle) = (DYSearly_diastolic_seg1(cycle) + DYSearly_diastolic_seg2(cycle) + DYSearly_diastolic_seg3(cycle) + DYSearly_diastolic_seg4(cycle) + DYSearly_diastolic_seg5(cycle))/5;
early_diastole_DYS = mean(DYS_DIAearly);
DYS_DIAlate (cycle) = (DYSlate_diastolic_seg1(cycle) + DYSlate_diastolic_seg2(cycle) + DYSlate_diastolic_seg3(cycle) + DYSlate_diastolic_seg4(cycle) + DYSlate_diastolic_seg5(cycle))/5;
late_diastole_DYS = mean(DYS_DIAlate);

% Calculate the entropy in the ECG signal
ENTROPY_ECG (cycle) = entropy (ecg_cycles);


%% SDI
cycle_duration = length(volume_cycles)*dt;
[~,t_minLVvolidx] = min(volume_cycles);
t_minLVvol = t_minLVvolidx*dt;
[~,t_minSeg1idx] = min(segm1_cycles);
t_minseg1 = t_minSeg1idx*dt;
[~,t_minSeg2idx] = min(segm2_cycles);
t_minseg2 = t_minSeg2idx*dt;
[~,t_minSeg3idx] = min(segm3_cycles);
t_minseg3 = t_minSeg3idx*dt;
[~,t_minSeg4idx] = min(segm4_cycles);
t_minseg4 = t_minSeg4idx*dt;
[~,t_minSeg5idx] = min(segm5_cycles);
t_minseg5 = t_minSeg5idx*dt;
t_seg1 = (t_minseg1/cycle_duration);
t_seg2 = (t_minseg2/cycle_duration);
t_seg3 = (t_minseg3/cycle_duration);
t_seg4 = (t_minseg4/cycle_duration);
t_seg5 = (t_minseg5/cycle_duration);

% differences for the new SDI
dif_segm1 = abs(t_minseg1-t_minLVvol);
dif_segm2 = abs(t_minseg2-t_minLVvol);
dif_segm3 = abs(t_minseg3-t_minLVvol);
dif_segm4 = abs(t_minseg4-t_minLVvol);
dif_segm5 = abs(t_minseg5-t_minLVvol);

SDI_cycle (cycle) = std([t_seg1 t_seg2 t_seg3 t_seg4 t_seg5]);
SDI  = mean(SDI_cycle)*100;
SDInew_cycle (cycle) = std([dif_segm1 dif_segm2 dif_segm3 dif_segm4 dif_segm5]);
SDInew = mean(SDInew_cycle)*100;

figure (42)
    plot(IF_cycles,'linewidth',0.5)
    hold on 
end

%% QRS analysis based on MODWT (Matlab documentation)
% Resample signal to 1Khz
targetSampleRate = 1e+03;
[sample_ECG,time_ECG] = resample(ECG_INCA,tt1_interp,targetSampleRate);
wt = modwt(sample_ECG,5);
wtrec = zeros(size(wt));
wtrec(4:5,:) = wt(4:5,:);
y = imodwt(wtrec,'sym4');
y = abs(y).^2;
thresholdx = max(y)*0.005;
[qrspeaks,locs] = findpeaks(y,time_ECG,'MinPeakHeight',thresholdx);
min_distance = length(qrspeaks)/length(Rwaves);
[rpeaks,locs_Rpeaks] = findpeaks(qrspeaks,'MinPeakDistance',min_distance*0.65);
R_waves = locs(locs_Rpeaks);
if length(locs(locs>R_waves(end))) < 2
    locs_Rpeaks = locs_Rpeaks(1:end-1);
    R_waves = R_waves(1:end-1);
end
if length(locs(locs<R_waves(1))) < 2
    locs_Rpeaks = locs_Rpeaks(2:end);
    R_waves = R_waves(2:end);
end
Q_waves = locs(locs_Rpeaks-2);
S_waves = locs(locs_Rpeaks+2);

QRS_durations_cycle = S_waves - Q_waves;
QRS_duration = mean(rmoutliers(QRS_durations_cycle,'median'));

figure (1001)
plot(time_ECG,y,'DisplayName','y');hold on;plot(locs,qrspeaks,'*');hold on; yyaxis right;plot(time_ECG,sample_ECG,'DisplayName','ecg_tonom');hold off

SDI_calculation;