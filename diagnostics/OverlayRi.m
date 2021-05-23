% Script to overlay Ri on wave travel diagram

% calculate BV freq (change code if C is changing with time)
BVf = BV(g(:,1),C(:,1),gamma(:,1));

for i = 1:length(T_arr)
Ri(:,i) = RichardsonNumber(BVf,LastDomainZindex,U(:,length(X_KM)/2,i) - wind(3:LastDomainZindex,length(X_KM)/2),dz);
end

% plot wave travel diagram
figure % wave travel plot for x in the middle of domain
contourf(T_arr,Z_KM,squeeze(W(:,length(X_KM)/2,:)).*SCALING_FACTOR(:,length(X_KM/2)),50,'Edgecolor','none')
xlabel('time (s)','Interpreter','latex','FontSize',12)
ylabel('z (km)','Interpreter','latex','FontSize',12)
hold on
colorbar
% overlay contours where Ri <= 0.25
contour(T_arr,Z_KM,Ri<=0.25,[0 1])
cmocean gray

%%save as pdf using MATLAB  exportgraphics function:
%exportgraphics(gcf,'./WTD_breaking.pdf','ContentType','image') 
%%change 'image' to 'vector' if needed -> a vector will give large size pdf
%%If vector image is too large and matlab can't generate a nicer looking image, use export_fig tool 
% export_fig FILENAME.png -png -r300 -transparent -painters % a 300 dpi PNG