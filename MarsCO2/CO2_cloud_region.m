
Tcond = 90 + ((104.5-Z_KM)./1.853);

Tcond = repmat(Tcond',[1,length(X_KM)]);

% overall temperature 
T = T0(3:LastDomainZindex,3:end-2)+ T_PERT;

% adjust z to reflect topography
Z_abs = Z + hq.*1000;

figure
% for pass 390
% contourf(X(3:LastDomainZindex,3:end-2)./1000,Z_abs(3:LastDomainZindex,3:end-2)./1000,Tcond-squeeze(T(:,:,end)),50,'Edgecolor','none')
% hold on
% plot(X_KM,hq(3:end-2)+0.28,'LineWidth',2,'color','black') % plot topography
% load('Pass390_cloudecho.mat')
% plot(Pass390_cloudecho(:,1),Pass390_cloudecho(:,2)+0.1,'.','MarkerSize',8,'color','black') % plot cloud echoes 


% For Pass 260
% contourf(X(3:LastDomainZindex,3:end-2)./1000,Z_abs(3:LastDomainZindex,3:end-2)./1000,Tcond-squeeze(T(:,:,end)),50,'Edgecolor','none')
% hold on
% plot(X_KM,hq(3:end-2)+0.28,'LineWidth',2,'color','black') % plot topography
% load('Pass260_cloudecho.mat') % load cloud echo data points
% Pass260_cloudecho(:,2) = Pass260_cloudecho(:,2)+0.5;
% plot(Pass260_cloudecho(:,1),Pass260_cloudecho(:,2),'.','MarkerSize',8,'color','black') % plot cloud echoes 

% For Pass 207
contourf(X(3:LastDomainZindex,3:end-2)./1000,Z_abs(3:LastDomainZindex,3:end-2)./1000,Tcond-squeeze(T(:,:,end)),50,'Edgecolor','none')
hold on
plot(X_KM,hq(3:end-2)+0.1,'LineWidth',2,'color','black') % plot topography
load('Pass207_echo.mat') % load cloud echo data points
Pass207_echo(Pass207_echo(:,2) < -2.1,:) = [];
plot(Pass207_echo(:,1),Pass207_echo(:,2),'.','MarkerSize',8,'color','black') % plot cloud echoes 
cmocean thermal
colorbar
ylim([-5 10])
xlim([0 550])

% for separate contour and terrain plot:
% figure
% subplot(2,1,1)
% contourf(X_KM,Z_KM,Tcond-squeeze(T(:,:,end)),50,'Edgecolor','none')
% hold on
% load('Pass390_cloudecho.mat')
% plot(Pass390_cloudecho(:,1),Pass390_cloudecho(:,2)+1,'+','color','black')
% ylim([0.25 3])
% colorbar
% caxis([-1 2])
% 
% subplot(2,1,2)
% plot(X_KM,hq(3:end-2),'LineWidth',1,'color','black')

% For ref: To find distance in km between 2 points on Mars: deg2km(distance(0,1,1,1),'mars')