clear

% forcing parameters
w0 = 0.001; % doesn't matter
omega = 0.007; % rad/s
T = 2*pi/omega;
t0 = 2400; % sec
sig_t = 600; % sec

dt = 1;
% time array to compute forcing function
t = 0:dt:12000; % sec

% forcing function (at x = 0, z = 0)
w = w0.*cos(omega.*(t-t0)).*exp(-((t-t0).^2)./(2*sig_t^2));

% ---------------------------

W = fft(w);
n = length(t);
fs = 1/dt;
f = (0:n-1)*(fs/n);     % frequency range (in Hz)
power = abs(W).^2/n;    % power of the DFT

figure()
subplot(1,2,1)
plot(t./60,w,'-k','LineWidth',1)
xlabel('$t$ (min)','Interpreter','latex','FontSize',12)
ylabel('$w$ (m/s)','Interpreter','latex','FontSize',12)
title('Forcing','Interpreter','latex','FontSize',14)

subplot(1,2,2)
plot(2*pi*f,power,'color',[0 0.5 0],'Linewidth',2) %(plotting f in rad/s)
xlim([0 0.02])
xlabel('$\omega$ (rad/s)','Interpreter','latex','FontSize',12)
ylabel('$|P(f(t))|^2$','Interpreter','latex','FontSize',12)
title('Frequency spectrum','Interpreter','latex','FontSize',14)