dt = 60;
t_end = 7200;
x = -5000000:10000:5000000; % m

% -- Tsunami parameters --
tau = 600; %time taken to grow to max amplitude (s)
c = -200; % Tsunami phase speed (m/s)
%sig = 20000; % parameter to describe Lx (m)
X = x./100000; % units of 100 km 
h0 = 1; % m
w0 = 9.5e-6; % coefficient to get desired value of w
u0 = 0; % wind speed at sea surface

%X = X-12;
for t = 0:dt:t_end
    
        if t < tau
    
        h = (t/tau).*(h0.*airy(1-X).*(X./2).*exp((2-X)./2)); % ref Wu et al, Frit et al
        w = (t/tau).*(w0).*(u0-c).*(h0/2).*exp((2-X)./2).*(airy(1-X)-airy(1,1-X).*X-airy(1-X).*(X./2));
        
        else
            X_steady = X-(c*(t-tau))/100000; % in 100 km units
            h = h0.*airy(1-X_steady).*(X_steady./2).*exp((2-X_steady)./2);
            w = (w0).*(u0-c).*(h0/2).*exp((2-X_steady)./2).*(airy(1-X_steady)-airy(1,1-X_steady).*X-airy(1-X_steady).*(X_steady./2));
        end
   

end

% To plot: set t = tau, then compute h ad w -> run the plot script below
    
    figure
    plot(x./1000,w.*1000,'Linewidth',1)
    ylabel('Vertical velocity (mm/s)','Interpreter','latex','FontSize',12)
    ylim([-1.2 1.2])
    
    yyaxis right
    plot(x./1000,h,'--','Linewidth',1)
    ylim([-1 1])
    ylabel('sea surface height (m)','Interpreter','latex','FontSize',12)
    xlabel('x (km)','Interpreter','latex','FontSize',12)
    xlim([-2500 2500])
    title('Tsunami forcing','Interpreter','latex','FontSize',14)
    %pause(0.1)
    






