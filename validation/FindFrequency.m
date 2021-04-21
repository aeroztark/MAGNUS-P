%This function will plot frequency spectrogram using Wavelet transform for a time series  

function [] = FindFrequency(t,y,dt)

%INPUTS
% t: time value of timeseries data in sec
% y: y value of timeseries data
% dt: data cadence in sec (15s in GPS TEC study)

%Functions used: wavelet.m and wave_signif.m : read docs for optional inputs

%If reading the timeseries from a .fig file, follow these lines of code
% fig = gcf;  %select current figure
% axsObjs = fig.Children; %get children
% dataObjs = axsObjs.Children;
% x = dataObjs(1).XData;  %assign data
% y = dataObjs(1).YData;

% Wavelet computation
[wave,period,scale,coi] = wavelet(y,dt);
power = (abs(wave)).^2 ;        % compute wavelet power spectrum

% Significance levels: (if computing)
% [signif,fft_theor] = wave_signif(y,dt,scale);
% sig95 = (signif')*(ones(1,length(y)));  % expand signif --> (J+1)x(N) array
% sig95 = power ./ sig95;         % where ratio > 1, power is significant

%Plotting

figure
contourf(t./60,1000./period,power,50,'EdgeColor','none');  %1000/period to plot freq on y axis in mHz
colorbar
ylabel('frequency (mHz)')
xlabel('time (min)')
cmocean thermal
%plot 95% confidence level (significance at 5%) contour
% hold on
% contour(x,1000./period,sig95,[-99,1],'k','LineWidth',1)

%plot COI
%plot(x,1000./coi,'k','LineWidth',2) %freq converted to mHz

