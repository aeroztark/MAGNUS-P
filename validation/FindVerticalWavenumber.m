function[Z_LOCATION, MAX_Lz] = FindVerticalWavenumber(profile, Z_KM, dz)

    %%% This function will take in a vertical profile (e.g. vertical velocity)
    %%% and will return the dominant vertical wavelength. This can be used to
    %%% validate the GW numerical solutions with predictions from linear dispersion relation. 
    
    %%% profile: 1D array of vertical profile
    %%% dz: resolution along vertical axis (in km)
    %%% Z_KM: array of vertical coordinate values in km
    %%% OUT: Z_LOCATION -> where max Lz is found; MAX_Lz -> value of max Lz
    
    [wave,period,scale,coi] = wavelet(profile,dz);  % Terrence & Compo's function
    power = (abs(wave)).^2 ;        % compute wavelet power spectrum (units sigma^2 (variance) of dz)
    
    % find max wavelet power
    M = max(power,[],1,'omitnan');
    [max_col_val,max_col_ind] = max(M); %find the column with max power
    [max_row_val,max_row_ind] = max(power(:,max_col_ind)); 
    
    Z_LOCATION = Z_KM(max_col_ind);
    MAX_Lz = period(max_row_ind);
    
    % Generate plot
    figure
    contourf(period,Z_KM,power',50,'EdgeColor','none') 
    xlabel('Vertical Wavelength (km)')
    ylabel('z (km)')
    title('Wavelet spectrum for vertical scale')
    colormap jet