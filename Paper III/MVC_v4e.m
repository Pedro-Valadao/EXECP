%% 
% Script for MVC EMG RMS analysis
% 
% F Cenni 6 Apr 2020, francesco.cenni@gmail.com. 
%
% Description: The script finds the maximum RMS EMG for Soleus and Medial Gastrocnemius 
% during the isometric part of the muscle action. Three sliding window lengths 
% were used: 50,100 and 200 ms.
%
% Last updated 3 Mar 2021
% correction on the method for calculating rms using movmean function,
% following this example https://it.mathworks.com/matlabcentral/answers/310418-how-do-i-calculate-the-rms-of-every-5-values
% 9 Mar 2021 -> movmean Return Only Full-Window Averages, M =
% movmean(A,3,'Endpoints','discard')... check movmean for peak_emgrms_sol_200_new2 (as example)
%pv:changed figure 4, fourth plot to MG

%%
clear all, close all, clc
dataPath = 'C:\Doutorado\Matlab\Paper II\MVC\smr_MVC\';
fileName = '5208_pre1_MVC.smr'; 
data = ImportSMR([dataPath fileName]);
Subject = string(fileName(1:4));
Test = "Pre2"; 
Group = "CP";
%%
%%
%Data Trimming 
data(1).imp.adc = data(1).imp.adc(121500:end);
data(2).imp.adc = data(2).imp.adc(121500:end);
data(3).imp.adc = data(3).imp.adc(2*121500:end);
data(4).imp.adc = data(4).imp.adc(2*121500:end);
%data(7).imp.adc = data(7).imp.adc(1:135000);
%data(8).imp.adc = data(8).imp.adc(1:135000);

%% channel_1:Torque data [Nm]
% information extracted from the struct
f_sample_ch1=1/(data(1).hdr.adc.SampleInterval(1)*data(1).hdr.adc.SampleInterval(2));
variable_name_ch1=data(1).hdr.title;  %name of the variable
scale_factor_ch1=data(1).hdr.adc.Scale; %scaling factor for channel1
unit_ch1=data(1).hdr.adc.Units; %unit channel1
channel1=data(1).imp.adc; % signal (format:int16)
N_samples_ch1=length(channel1); %number of samples
channel1_scaled=(double(channel1))*scale_factor_ch1; %scaling
time = N_samples_ch1/f_sample_ch1; %s
fprintf('the measurements took %d s\n',time); % based on the TRUSTED freq of channel1

%filter the data in channel1 (torque)
[b,a]=butter(6,40/(f_sample_ch1/2),'low'); %low-pass filter 40Hz
%channel1_scaled_d = double(channel1_scaled); %conversion to double
channel1_scaled_filt = filtfilt(b,a,channel1_scaled); %signal filtered
channel1_scaled_window = movmean(channel1_scaled,20); %where each mean is calculated over a sliding window 
% offset detection
starting_value_torque=mean(channel1_scaled_window([1:200])); %average from first 200 samples (arbitrarily chosen)
% offset removed for the torque data
torque_data = channel1_scaled_window - starting_value_torque;

%% channel_2 (Angle)
f_sample_ch2 = 1/(data(2).hdr.adc.SampleInterval(1)*data(2).hdr.adc.SampleInterval(2));
variable_name_ch2=data(2).hdr.title;  
scale_factor_ch2=data(2).hdr.adc.Scale;
unit_ch2=data(2).hdr.adc.Units; %unit channel2
channel2=data(2).imp.adc;  
N_samples_ch2=length(channel2);
channel2_scaled_int=channel2*scale_factor_ch2; 
channel2_scaled_d=(double(channel2))*scale_factor_ch2;% double format 

%% pre-plots
figure(1)
plot(channel2_scaled_int)
ylabel('angles ')
margin_a = input('Enter interval for calculating starting angle before stretch, if no artefacts select 1-200, otherwise here inferior margin in frame: ','s');
margin_b = input('here superior margin in frame: ','s');
margin_a = str2double(margin_a);
margin_b = str2double(margin_b);

% offset detection
starting_value_int=int8(mean(channel2_scaled_int([margin_a:margin_b])));
%starting_value_int=mean(channel2_scaled_int([100000:100200])); %%---> only for 5117
plantarflex_value_int=starting_value_int+28; %prior knowledge: always 30 deg in plantarflexion???
%starting_value_d=mean(channel2_scaled_d([1:200])); %average from first 200 samples
starting_value_d=mean(channel2_scaled_d([margin_a:margin_b])); %%---> only for 5117
plantarflex_value_d=starting_value_d+28; %

%pedal pre-signal processing (needed because there are some spikes on the angle signal +/-1 deg)
%indices
i=1;j=1;l=1;m=1;
for c = 1:N_samples_ch2
    if c == N_samples_ch2-200 % to avoid to exceed the N_samples_ch2 array
        break
    end
    if channel2_scaled_int(c) > plantarflex_value_int;
        channel2_scaled_int(c)=plantarflex_value_int;
    end
    if channel2_scaled_int(c) < starting_value_int;
        channel2_scaled_int(c)=starting_value_int;
    end
    if channel2_scaled_int(c) == starting_value_int+1;
        if channel2_scaled_int(c+200) < starting_value_int+2;
            channel2_scaled_int(c)=starting_value_int;
        end
    end
end
channel2_scaled_int = int8(medfilt1(double(channel2_scaled_int),2000)); %further clean the signal
%-->exception for 5119
%channel2_scaled_int = int8(medfilt1(double(channel2_scaled_int(30000:end,2)),2000)); %further clean the signal

%starting_value_int=int8(mean(channel2_scaled_int([margin_a:margin_b])));
max_plantar_values_int= max(channel2_scaled_int);
%plantarflex_value=starting_value+30; %prior knowledge 20deg plantar
plantarflex_value_int = int8(max_plantar_values_int); %from real max value

%-->DETECTING the EVENTS from the angle signal
for c = 1:N_samples_ch2 %exception for 5119 add--> -30000
    if channel2_scaled_int(c) == starting_value_int+1;
        if channel2_scaled_int(c-1) == starting_value_int
            %first event trapezoidal wave form (ankle is started to be plantarflexed)    
            index_ev1(i)=c-1;
            i=i+1;
        end
        if channel2_scaled_int(c+1) == starting_value_int
            %forth event trapez wave form (ankle is now in dorsiflexion)    
            index_ev4(j)=c+1;
            j=j+1;
        end
    end
    if channel2_scaled_int(c) == plantarflex_value_int;
        if channel2_scaled_int(c-1) == plantarflex_value_int-1;
            %second event trapez wave form (ankle is now in plantarflexion)    
            index_ev2(l)=c;
            l=l+1;
        end
        if channel2_scaled_int(c+1) == plantarflex_value_int-1;
            %third event trapez wave form (ankle is started to be dorsiflexed)  
            index_ev3(m)=c;
            m=m+1;
        end
    end
end

% offset for the angle data (channel2_scaled_shifted)
angle_data_int = channel2_scaled_int - starting_value_int;
angle_data_d = channel2_scaled_d - starting_value_d;

%% EMG data

%% channels avalable for wired system
channel3=data(3).imp.adc; %soleus
channel4=data(4).imp.adc; %MG
channel3_d=double(channel3);%
channel4_d=double(channel4);%
%% channels avalable for wireless system
channel7=data(7).imp.adc; %soleus
channel8=data(8).imp.adc; %MG
channel7_d=double(channel7);%
channel8_d=double(channel8);%

%low-pass filter (freq=4KHz)
fcuthigh=195; %high cut frequency in Hz
[b1,a1]=butter(6,fcuthigh/(4000/2),'low'); %order=6, %fs=4000 
%bandpass
fcutlow=20;   %low cut frequency in Hz
fcuthigh=195;   %high cut frequency in Hz
[b0,a0]=butter(6,[fcutlow,fcuthigh]/(4000/2),'bandpass'); %order=6, %fs=4000

%% pre-plots
figure(1)
subplot(2,1,1)
plot(channel3)
ylabel('Soleus EMG ')
title('wired system')
hold on
subplot(2,1,2)
plot(channel4,'r')
ylabel('Medial Gastrocnemius EMG ')

figure(2)
subplot(2,1,1)
plot(channel7)
ylabel('Soleus EMG ')
title('wireless system')
hold on
subplot(2,1,2)
plot(channel8,'r')
ylabel('Medial Gastrocnemius EMG ')
%% condition for selecting the type of EMG channel

typeEMG = input('Enter type EMG: ','s');

switch typeEMG
    case 'wired'
        %------channel 3 (WIRED channel SOLEUS)
        variable_name_ch3=data(3).hdr.title;  
        scale_factor_ch3=data(3).hdr.adc.Scale; 
        unit_ch3=data(3).hdr.adc.Units; %unit channel3
        channel3=data(3).imp.adc; % signal as int16
        N_samples_ch3=length(channel3);
        channel3_d=double(channel3);
        channel3_scaled=channel3_d*scale_factor_ch3; %scaled emg as double
        freq_ratio = int8(N_samples_ch3./N_samples_ch2);  
        %resampling of the emg signal with respect the one from the angle 
        channel3_d2 = resample(channel3_scaled,1,double(freq_ratio)); 
        %Use -filtfilt- when it is critical to keep the PHASE information of a signal INTACT
        channel3_filt = filtfilt(b1,a1,channel3_d2); %signal filtered LOW PASS
        channel3_filt2 = filtfilt(b0,a0,channel3_d2); %signal filtered BAND PASS

        %----channel 4 (WIRED channel MEDIAL GASTROCNEMIUS)
        variable_name_ch4=data(4).hdr.title;  
        scale_factor_ch4=data(4).hdr.adc.Scale; %
        unit_ch4=data(4).hdr.adc.Units; %unit channel4
        channel4=data(4).imp.adc; 
        N_samples_ch4=length(channel4);
        channel4_d=double(channel4);
        channel4_scaled=channel4_d*scale_factor_ch4; %scaled emg as double
        %resampling of the emg signal with respect the one from the angle 
        channel4_d2 = resample(channel4_scaled,1,double(freq_ratio)); %channel4_d
        %Use -filtfilt- when it is critical to keep the PHASE information of a signal INTACT
        channel4_filt = filtfilt(b1,a1,channel4_d2); %signal filtered LOW PASS
        channel4_filt2 = filtfilt(b0,a0,channel4_d2); %signal filtered BAND PASS
        
    case 'wireless'
        %------channel 7 (wireless channel SOLEUS) 5101 = ch 6!!
        t_delay_ch7 = 200; %ms !!!-->the emg data have a constant delay of 200ms with respect angle data
        variable_name_ch3=data(7).hdr.title;  
        scale_factor_ch3=data(7).hdr.adc.Scale; %this seems not needed
        unit_ch3=data(7).hdr.adc.Units; %unit channel7
        channel3=data(7).imp.adc; % signal as int16
        N_samples_ch3=length(channel3);
        channel3_d=double(channel3);
        channel3_scaled=channel3_d*scale_factor_ch3;  %scaled emg as double
        freq_ratio = int8(N_samples_ch3./N_samples_ch2);  
        %resampling of the emg signal with respect the one from the angle 
        channel3_d2 = resample(channel3_scaled,1,double(freq_ratio));
        %Use -filtfilt- when it is critical to keep the PHASE information of a signal INTACT
        channel3_filt = filtfilt(b1,a1,channel3_d2); %signal filtered LOW PASS
        channel3_filt2 = filtfilt(b0,a0,channel3_d2); %signal filtered BAND PASS

        %----channel 8 (wireless channel GASTROCNEMIUS)
        t_delay_ch8 = 200; %ms !!!-->the emg data have a constant delay of 200ms with respect angle data
        variable_name_ch4=data(8).hdr.title;  
        scale_factor_ch4=data(8).hdr.adc.Scale; %
        unit_ch4=data(8).hdr.adc.Units; %unit channel8
        channel4=data(8).imp.adc; 
        N_samples_ch4=length(channel4);
        channel4_d=double(channel4);
        channel4_scaled=channel4_d*scale_factor_ch4; %scaled emg as double
        %resampling of the emg signal with respect the one from the angle 
        channel4_d2 = resample(channel4_scaled,1,double(freq_ratio)); %channel8_d
        %Use -filtfilt- when it is critical to keep the PHASE information of a signal INTACT
        channel4_filt = filtfilt(b1,a1,channel4_d2); %signal filtered LOW PASS
        channel4_filt2 = filtfilt(b0,a0,channel4_d2); %signal filtered BAND PASS        

end        
close all

%% select relevant events (conditions: onset movement 'ev1' and torque threshold ...>2Nm)

cc=1;dd=1;
threshold_torque = 2;
t_sel =15; %search intial torque value 15s before onset movement
t_back =5; %time in s for background noise
for ii = 1:length(index_ev1) %
        if torque_data(index_ev1(ii)) > threshold_torque 
            events_torque_sol_sel(cc)=index_ev1(ii);
            cc = cc+1;
            for jj = (index_ev1(ii) - t_sel*f_sample_ch1):10:index_ev1(ii) %every 10 samples
                if torque_data(jj) > 10
                    events_init_torque_sol_sel(dd)=jj;
                    dd=dd+1;
                    break
                end
            end
            if dd < cc
                events_init_torque_sol_sel(dd)=0;
                dd=dd+1;
            end
        end            
            
end

%% peak rms calculations
%SOL
% using this approach, the mean for rms is outside sqrt
emgrms_sol = sqrt(abs(channel3_filt2.^2)+1E-14); %rectified emg SOL (rms just to keep old name)
emgrms_sol_mw_50 = movmean(emgrms_sol,50); %%where each mean is calculated over a sliding window 
%CORRECTION 2/3/21: bandpass signal --> squared, avg based on win size and
%sqrt; this is performed over the filtered bandpass signal
emgrms_sol_mw_200_new = sqrt(movmean(channel3_filt2.^2,200)); %

%MG
emgrms_mg = sqrt(abs(channel4_filt2.^2)+1E-14); % rectified emg MG (rms just to keep old name)
emgrms_mg_mw_200 = movmean(emgrms_mg,200); %%where each mean is calculated over a sliding window 
%CORRECTION 2/3/21: bandpass signal --> squared, avg based on win size and
%sqrt; this is performed over the entire (filtered bandpass) signal
emgrms_mg_mw_200_new = sqrt(movmean(channel4_filt2.^2,200));

ic=1;
for ii = 1:length(events_init_torque_sol_sel)
    if events_init_torque_sol_sel(ii) ~= 0
        events_init_torque_sol_sel_n(ic) =  events_init_torque_sol_sel(ii);
        events_torque_sol_sel_n(ic) = events_torque_sol_sel(ii);
        ic=ic+1;
    end
end
        

for ii = 1:length(events_init_torque_sol_sel_n)
    switch typeEMG
        case 'wired'
        %if typeEMG == 'wired'
            %SOL        
            %peak_emgrms_sol_200(ii) = max(emgrms_sol_mw_200(events_init_torque_sol_sel(ii):events_torque_sol_sel(ii))); %old
            %peak_emgrms_sol_200_new(ii) = max(emgrms_sol_mw_200_new(events_init_torque_sol_sel(ii):events_torque_sol_sel(ii)));
            % sol signal extracted for the relevant section
            emg_sol_cut = channel3_filt2(events_init_torque_sol_sel_n(ii):(events_torque_sol_sel_n(ii)));
            % M = movmean(A,3,'Endpoints','discard') ONLY full window average
            peak_emgrms_sol_200_new2(ii) = max(sqrt(movmean(emg_sol_cut.^2,200,'Endpoints','discard')));
            %MG
            % mg signal extracted for the relevant section
            emg_mg_cut = channel4_filt2(events_init_torque_sol_sel_n(ii):(events_torque_sol_sel_n(ii)));
            % M = movmean(A,3,'Endpoints','discard') ONLY full window average
            peak_emgrms_mg_200_new2(ii) = max(sqrt(movmean(emg_mg_cut.^2,200,'Endpoints','discard')));
        %else
        case 'wireless'
        %SOL
            % sol signal extracted for the relevant section
            emg_sol_cut = channel3_filt2(events_init_torque_sol_sel_n(ii)-t_delay_ch7:(events_torque_sol_sel_n(ii)-t_delay_ch7));
            % M = movmean(A,3,'Endpoints','discard') ONLY full window average
            peak_emgrms_sol_200_new2(ii) = max(sqrt(movmean(emg_sol_cut.^2,200,'Endpoints','discard')));
            %MG
            % mg signal extracted for the relevant section
            emg_mg_cut = channel4_filt2(events_init_torque_sol_sel_n(ii)-t_delay_ch8:(events_torque_sol_sel_n(ii)-t_delay_ch8));
            % M = movmean(A,3,'Endpoints','discard') ONLY full window average
            peak_emgrms_mg_200_new2(ii) = max(sqrt(movmean(emg_mg_cut.^2,200,'Endpoints','discard')));
    end
end
%SOL     
fprintf('NEW maximum RMS-SOL valueS in a 200ms sliding window %.4f V\n',peak_emgrms_sol_200_new2);
max_peak_emgrms_sol_200_new = max(peak_emgrms_sol_200_new2);
fprintf('NEW max among trials of maximum RMS-SOL value in a 200ms sliding window %.4f V\n',max_peak_emgrms_sol_200_new);
%MG
fprintf('NEW maximum RMS-MG valueS in a 200ms sliding window %.4f V\n',peak_emgrms_mg_200_new2);
max_peak_emgrms_mg_200_new = max(peak_emgrms_mg_200_new2);
fprintf('NEW max among trials of maximum RMS-MG value in a 200ms sliding window %.4f V\n',max_peak_emgrms_mg_200_new);

% t_back --> last 5s before MVC
for ii = 1:length(index_ev1)
    switch typeEMG
        case 'wired'
            %SOL       
            p2p_filt_sol(ii) = peak2peak(channel3_filt2((index_ev1(ii)-t_back*f_sample_ch1):index_ev1(ii)));%on bandpass signal        
            emgrms_background_noise_sol=emgrms_sol((index_ev1(ii)-t_back*f_sample_ch1):index_ev1(ii)); %old
            avgs_emgrms_background_noise_sol(ii)=mean(emgrms_background_noise_sol); %old
            emgrms_background_noise_sol_new(ii)=rms(channel3_filt2((index_ev1(ii)-t_back*f_sample_ch1):index_ev1(ii))); %corrected
            %MG
            p2p_filt_mg(ii) = peak2peak(channel4_filt2((index_ev1(ii)-t_back*f_sample_ch1):index_ev1(ii)));
            emgrms_background_noise_mg=emgrms_mg((index_ev1(ii)-t_back*f_sample_ch1):index_ev1(ii)); %old
            avgs_emgrms_background_noise_mg(ii)=mean(emgrms_background_noise_mg); %old
            emgrms_background_noise_mg_new(ii)=rms(channel4_filt2((index_ev1(ii)-t_back*f_sample_ch1):index_ev1(ii))); %corrected
        case 'wireless'
            %SOL
            p2p_filt_sol(ii) = peak2peak(channel3_filt2((index_ev1(ii)-t_back*f_sample_ch1-t_delay_ch7):index_ev1(ii)-t_delay_ch7));   
            emgrms_background_noise_sol=emgrms_sol((index_ev1(ii)-t_back*f_sample_ch1-t_delay_ch7):index_ev1(ii)-t_delay_ch7); %last 5s before MVC
            avgs_emgrms_background_noise_sol(ii)=mean(emgrms_background_noise_sol);
            emgrms_background_noise_sol_new(ii)=rms(channel3_filt2((index_ev1(ii)-t_back*f_sample_ch1-t_delay_ch7):index_ev1(ii)-t_delay_ch7)); %last 5s before MVC
            %MG
            p2p_filt_mg(ii) = peak2peak(channel4_filt2((index_ev1(ii)-t_back*f_sample_ch1-t_delay_ch7):index_ev1(ii)-t_delay_ch7)); 
            emgrms_background_noise_mg=emgrms_mg((index_ev1(ii)-t_back*f_sample_ch1-t_delay_ch7):index_ev1(ii)-t_delay_ch7); %last 5s before MVC
            avgs_emgrms_background_noise_mg(ii)=mean(emgrms_background_noise_mg); 
            emgrms_background_noise_mg_new(ii)=rms(channel4_filt2((index_ev1(ii)-t_back*f_sample_ch1-t_delay_ch7):index_ev1(ii)-t_delay_ch7)); %last 5s before MVC

    end
    
    if p2p_filt_sol(ii) > 0.045 | emgrms_background_noise_sol_new > 0.05   %avgs_emgrms_background_noise_sol
        fprintf('SOL warning high background noise! trial%.2i \n',ii);
    end
    if p2p_filt_mg(ii) > 0.045 | emgrms_background_noise_mg_new > 0.05  %avgs_emgrms_background_noise_mg
        fprintf('MG warning high background noise! trial%.2i \n',ii);
    end
 
end
%% figures

figure(1)
subplot(2,1,1)
plot(channel1_scaled,'r')
ylabel('raw Ankle Torque [Nm]')
grid on
subplot(2,1,2)
hold on
plot(channel1_scaled_filt,'g') 
plot(channel1_scaled_window,'b')
legend({'filter','moving win'},'Location','northwest')
grid on
ylabel('Ankle Torque [Nm]')
xlabel('time [samples]')
% 
figure(2)
plot(channel3_d2,'r') 
hold on
%plot(channel3_filt,'g') %low pass
plot(channel3_filt2,'b') %band pass
grid on
legend({'raw','bandpass'},'Location','northwest')
ylabel('EMG SOL [mV]') 
xlabel('time [samples]')

figure(3)
plot(channel4_d2,'r') 
hold on
%plot(channel4_filt,'g') %low pass
plot(channel4_filt2,'b') %band pass
legend({'raw','bandpass'},'Location','northwest')
grid on
ylabel('EMG MG [mV]') 
xlabel('time [samples]')

figure(4)
subplot(4,1,1)
hold on
plot(angle_data_int) 
grid on
ylabel('angles [deg]') 
plot (index_ev1,0,'r.')
plot (index_ev2,plantarflex_value_int-starting_value_int,'k.')
plot (index_ev3,plantarflex_value_int-starting_value_int,'y.')
plot (index_ev4,0,'g.')
%
subplot(4,1,2)
hold on
plot(channel1_scaled_filt-starting_value_torque,'g') %signal filtered
plot(torque_data,'b') %signal obtained after moving win
plot(events_torque_sol_sel_n,torque_data(events_torque_sol_sel_n),'r.')
plot(events_init_torque_sol_sel_n,torque_data(events_init_torque_sol_sel_n),'r.')
grid on
ylabel('Ankle Torque [Nm]')
%
subplot(4,1,3)
hold on
plot(channel3_filt2) %channel3_d
plot(events_torque_sol_sel_n,0,'r.')
plot(events_init_torque_sol_sel_n,0,'r.')
%ylim([-0.3 0.3])
grid on
ylabel('EMG SOL [mV]')
%
subplot(4,1,4)
hold on
plot(channel4_filt2) 
plot(events_torque_sol_sel_n,0,'r.')
plot(events_init_torque_sol_sel_n,0,'r.')
%ylim([-0.3 0.3])
grid on
ylabel('EMG MG [mV]')

figure(5)
hold on
%plot(emgrms_sol,'b')
plot(emgrms_sol_mw_200_new,'k')
plot(events_torque_sol_sel_n,0,'r.')
plot(events_init_torque_sol_sel_n,0,'r.')
%plot(events_init_torque_sol_sel,peak_emgrms_sol_200,'ro')
plot(events_init_torque_sol_sel_n,peak_emgrms_sol_200_new2,'go')
%ylim([0 0.1])
grid on
ylabel('RMS EMG-SOL [V]')

figure(6)
hold on
%plot(emgrms_sol,'b')
plot(emgrms_mg_mw_200_new,'g')
plot(events_torque_sol_sel_n,0,'r.')
plot(events_init_torque_sol_sel_n,0,'r.')
%plot(events_init_torque_sol_sel+50,peak_emgrms_mg_200,'go')
plot(events_init_torque_sol_sel_n+50,peak_emgrms_mg_200_new2,'ro')
grid on
%ylim([0 0.1])
ylabel('RMS EMG-MG [V]')

figure(7)
hold on
plot(events_init_torque_sol_sel_n,peak_emgrms_sol_200_new2,'r*')
plot(events_init_torque_sol_sel_n+50,peak_emgrms_mg_200_new2,'g*')
plot(events_torque_sol_sel_n,0,'r.')
plot(events_init_torque_sol_sel_n,0,'r.')
grid on
ylabel('PEAK RMS EMG-SOL [V]')
legend({'SOL 200ms','MG 200ms'},'Location','southwest')

figure(8)
hold on
plot(emgrms_background_noise_mg_new,'b.')
plot(emgrms_background_noise_sol_new,'r.')
plot(p2p_filt_mg,'bo')
hold on
plot(p2p_filt_sol,'ro')
grid on
ylabel('peak-to-peak(o) and background noise RMS(.) EMG [V]')
%%
%Save output
MVC_results = table(Subject,Group,Test,max_peak_emgrms_sol_200_new, max_peak_emgrms_mg_200_new)
save(Subject, "MVC_results")
%Notes:
%5218 MG EMG is 10x amplified, so it must be divided by 10 (line 218)
%5120 SOL EMG gain is 10x less! Multiply result by 10. 