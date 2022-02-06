function [omega_radps, f_Hz, T_minutes] = DispersionRelation(Lx,Lz,N_rad/s)
    % Get omega from Lx,Lz (meters), N (rad/s) -> all scalars
    
    omega_radps = (N*(2*pi/Lx))/sqrt(((2*pi/Lx).^2 + (2*pi/Lz).^2));
    f_Hz = omega_radps./(2*pi);
    T_minutes = 1/(f_Hz*60);
    
end