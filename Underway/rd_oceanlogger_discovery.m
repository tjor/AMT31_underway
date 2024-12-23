function tmp = rd_oceanlogger_discovery(fn_surf, fn_met, fn_light, fn_tsg)

    # Function reads & combines meta data from 4 different files (surf, met, light, tsg) into a single stucture: tmp

    # This function was modfied by tjor (Aug 2022) - key changes:
    # (i) function takes 4 separate filenames as input arguement
    # (ii) internal loops removed (previously function acted on set of files for cruise). Note: loops just act on element 1,
    # but keeps previous snytax to avoid any errors.


    pkg load netcdf

    % 1) surf
    for inc = 1:1 % loop acts on element 1:
        ncfile = fn_surf; # temporary variable
        if inc == 1
            % Assumes time is in dats (matlab format)
            tmp.time = ncread(ncfile,'time') + datenum([1899,12,30,0,0,0]);
            tmp.flowrate = ncread(ncfile,'flow'); % Instrument Flow Rate [l/mn]
            tmp.fluo = ncread(ncfile,'fluo'); % Fluorescence [V]
            tmp.trans = ncread(ncfile,'trans'); % Transmissibility [V]
            % These are all 0; need to get them from TSG
            % tmp.thermosalinograph_temp = nc{'temp_h'}(:); % Housing temp [degC]
            % tmp.sst = nc{'temp_m'}(:); % Remote temperature [degC]
            % tmp.conductivity =nc{'cond'}(:); % [S/m]
        else
            _time = ncread(ncfile,'time')+datenum([1899,12,30,0,0,0]);
           % Add a NaN between the two series (so that interpolation leaves a gap)
           midtime = (_time(1)-tmp.time(end))/2;
           tmp.time = [tmp.time; tmp.time(end)+midtime; _time];
           tmp.flowrate = [tmp.flowrate; NaN; ncread(ncfile,'flow')];
           tmp.fluo = [tmp.fluo; NaN; ncread(ncfile,'fluo')];
           tmp.trans = [tmp.trans; NaN; ncread(ncfile,'trans')];
        endif
    endfor

    % 2) met
    for inc = 1:1
        ncfile = fn_met; # temporary variable
        if inc == 1
            % tmp2.time = nc{'time'}(:)+datenum([1899,12,30,0,0,0]);
            tmp.wind_vel = ncread(ncfile,'speed'); %WInd speed [m/s]
            tmp.air_temp = ncread(ncfile,'airtemp'); % Air temp [degC]
            tmp.wind_dir = ncread(ncfile,'direct'); % Wind dir [deg]
            tmp.humidity = ncread(ncfile,'humid'); % Rel air humidity [%]
        else
           tmp.wind_vel = [tmp.wind_vel; NaN; ncread(ncfile,'speed')];
           tmp.air_temp = [tmp.air_temp; NaN; ncread(ncfile,'airtemp')];
           tmp.wind_dir = [tmp.wind_dir; NaN; ncread(ncfile,'direct')];
           tmp.humidity = [tmp.humidity; NaN; ncread(ncfile,'humid')];
        endif
    endfor

    % 3) Light
    for inc = 1:1
        ncfile = fn_light; # temporary variable
        if inc == 1
            % tmp3.time = nc{'time'}(:)+datenum([1899,12,30,0,0,0]);
            tmp.baro = ncread(ncfile,'pres'); % Atmospheric pressure [mbar]
            tmp.par1 = ncread(ncfile,'ppar'); % Port PAR [volt x 10-4]
            tmp.tir1 = ncread(ncfile,'ptir'); % port total irradiance [volt x 10-4]
            tmp.par2 = ncread(ncfile,'spar'); % Starboard PAR [volt x 10-4]
            tmp.tir2 = ncread(ncfile,'stir'); % Starboard total irradiance [volt x 10-4]
        else
           tmp.baro = [tmp.baro; NaN; ncread(ncfile,'pres')];
           tmp.par1 = [tmp.par1; NaN; ncread(ncfile,'ppar')];
           tmp.tir1 = [tmp.tir1; NaN; ncread(ncfile,'ptir')];
           tmp.par2 = [tmp.par2; NaN; ncread(ncfile,'spar')];
           tmp.tir2 = [tmp.tir2; NaN; ncread(ncfile,'stir')];
        endif
    endfor

    % 4) TSG (different time => need interpolation)
    for inc = 1:1
        ncfile = fn_tsg;
        if inc == 1
            tmp2.time = ncread(ncfile,'time') + datenum([1899,12,30,0,0,0]);
            tmp2.sal = ncread(ncfile,'salin'); % TSG salinity
            tmp2.sst = ncread(ncfile,'temp_r'); % remote temp
            tmp2.thermosalinograph_temp = ncread(ncfile,'temp_h');
            tmp2.conductivity =ncread(ncfile,'cond'); % [S/m]
        else
            _time = ncread(ncfile,'time')+datenum([1899,12,30,0,0,0]);
           % Add a NaN between the two series (so that interpolation leaves a gap)
           midtime = (_time(1)-tmp2.time(end))/2;
           tmp2.time = [tmp2.time; tmp2.time(end)+midtime; _time];
           tmp2.sal = [tmp2.sal; NaN; ncread(ncfile,'salin')];
           tmp2.sst = [tmp2.sst; NaN; ncread(ncfile,'temp_r')];
           tmp2.thermosalinograph_temp = [tmp2.thermosalinograph_temp; NaN; ncread(ncfile,'temp_h')];
           tmp2.conductivity = [tmp2.conductivity; NaN; ncread(ncfile,'cond')];
        endif
    endfor

    % Interpolate TSG variables to SURF time
    tmp.sal = interp1(tmp2.time,tmp2.sal,tmp.time);
    tmp.sst = interp1(tmp2.time,tmp2.sst,tmp.time);
    tmp.thermosalinograph_temp = interp1(tmp2.time,tmp2.thermosalinograph_temp,tmp.time);
    tmp.conductivity = interp1(tmp2.time,tmp2.conductivity,tmp.time);

    %tmp.chl = d(:,24); % [ug/l]  ? tjor: optional variables? keep for now!
    %tmp.sample_temp = d(:,25); % [degC]
    %
    %tmp.field28th = d(:,28);
    %tmp.field29th = d(:,29);

endfunction
