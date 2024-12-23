function tmp = rd_oceanlogger_cook(fn_surf, fn_met, fn_light, fn_tsg)

    # Function reads & combines meta data from 4 different files (surf, met, light, tsg) into a single stucture: tmp

    # This function was modfied by tjor (Aug 2022) - key changes:
    # (i) function takes 4 separate filenames as input arguement
    # (ii) internal loops removed (previously function acted on set of files for cruise). Note: loops just act on element 1,
    # but keeps previous snytax to avoid any errors.


    # This function was also tested on James Cook Nov 2024. Mods:

    #  (i) Flow variable is `flow' on JC, rather than flow2 on Disco.
    #  (ii) Loops over ship metadata files were re-implemented (as often there were multiple files)

    pkg load netcdf

    % 1) surf
    for inc = 1:length(fn_surf) % loop acts on element 1:
        ncfile = fn_surf(inc); # temporary variable
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
    for inc = 1:length(fn_met)
        ncfile = fn_met(inc); # temporary variable
        if inc == 1
            tmp2.time = ncread(ncfile,'time') + datenum([1899,12,30,0,0,0]);
            tmp2.wind_vel = ncread(ncfile,'speed'); %WInd speed [m/s]
            tmp2.air_temp = ncread(ncfile,'airtemp'); % Air temp [degC]
            tmp2.wind_dir = ncread(ncfile,'direct'); % Wind dir [deg]
            tmp2.humidity = ncread(ncfile,'humid'); % Rel air humidity [%]
        else
           _time = ncread(ncfile,'time')+datenum([1899,12,30,0,0,0]);
           % Add a NaN between the two series (so that interpolation leaves a gap)
           midtime = (_time(1)-tmp2.time(end))/2;
           tmp2.time = [tmp2.time; tmp2.time(end)+midtime; _time];
           tmp2.wind_vel = [tmp2.wind_vel; NaN; ncread(ncfile,'speed')];
           tmp2.air_temp = [tmp2.air_temp; NaN; ncread(ncfile,'airtemp')];
           tmp2.wind_dir = [tmp2.wind_dir; NaN; ncread(ncfile,'direct')];
           tmp2.humidity = [tmp2.humidity; NaN; ncread(ncfile,'humid')];
        endif
    endfor

      % Interpolate met variables to SURF time
    tmp.wind_vel = interp1(tmp2.time,tmp2.wind_vel,tmp.time,'extrap');
    tmp.air_temp = interp1(tmp2.time,tmp2.air_temp,tmp.time,'extrap');
    tmp.wind_dir = interp1(tmp2.time,tmp2.wind_dir,tmp.time,'extrap');
    tmp.humidity = interp1(tmp2.time,tmp2.humidity,tmp.time,'extrap');

    tmp2 =[];

    % 3) Light
    for inc = 1:length(fn_light)
        ncfile = fn_light(inc); # temporary variable
        if inc == 1
            tmp2.time = ncread(ncfile,'time') + datenum([1899,12,30,0,0,0]);
            tmp2.baro = ncread(ncfile,'pres'); % Atmospheric pressure [mbar]
            tmp2.par1 = ncread(ncfile,'ppar'); % Port PAR [volt x 10-4]
            tmp2.tir1 = ncread(ncfile,'ptir'); % port total irradiance [volt x 10-4]
            tmp2.par2 = ncread(ncfile,'spar'); % Starboard PAR [volt x 10-4]
            tmp2.tir2 = ncread(ncfile,'stir'); % Starboard total irradiance [volt x 10-4]
        else
           _time = ncread(ncfile,'time')+datenum([1899,12,30,0,0,0]);
           % Add a NaN between the two series (so that interpolation leaves a gap)
           midtime = (_time(1)-tmp2.time(end))/2;
           tmp2.time = [tmp2.time; tmp2.time(end)+midtime; _time];
           tmp2.baro = [tmp2.baro; NaN; ncread(ncfile,'pres')];
           tmp2.par1 = [tmp2.par1; NaN; ncread(ncfile,'ppar')];
           tmp2.tir1 = [tmp2.tir1; NaN; ncread(ncfile,'ptir')];
           tmp2.par2 = [tmp2.par2; NaN; ncread(ncfile,'spar')];
           tmp2.tir2 = [tmp2.tir2; NaN; ncread(ncfile,'stir')];
        endif
    endfor

    tmp.baro = interp1(tmp2.time,tmp2.baro,tmp.time,'extrap');
    tmp.par1 = interp1(tmp2.time,tmp2.par1,tmp.time,'extrap');
    tmp.tir1 = interp1(tmp2.time,tmp2.tir1,tmp.time,'extrap');
    tmp.par2 = interp1(tmp2.time,tmp2.par2,tmp.time,'extrap');
    tmp.tir2 = interp1(tmp2.time,tmp2.tir2,tmp.time,'extrap');
    tmp2 =[];

    % 4) TSG (different time => need interpolation)
    for inc = 1:length(fn_tsg)
        ncfile = fn_tsg(inc);
        if inc == 1
            tmp2.time = ncread(ncfile,'time') + datenum([1899,12,30,0,0,0]);
            tmp2.sal = ncread(ncfile,'salin'); % TSG salinity
            tmp2.sst = ncread(ncfile,'temp_r'); % remote temp
            tmp2.thermosalinograph_temp = ncread(ncfile,'temp_h');
            tmp2.conductivity = ncread(ncfile,'cond'); % [S/m]
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
    tmp.sal = interp1(tmp2.time,tmp2.sal,tmp.time,'extrap');
    tmp.sst = interp1(tmp2.time,tmp2.sst,tmp.time,'extrap');
    tmp.thermosalinograph_temp = interp1(tmp2.time,tmp2.thermosalinograph_temp,tmp.time,'extrap');
    tmp.conductivity = interp1(tmp2.time,tmp2.conductivity,tmp.time,'extrap');

   # keyboard
    %tmp.chl = d(:,24); % [ug/l]  ? tjor: optional variables? keep for now!
    %tmp.sample_temp = d(:,25); % [degC]

    %tmp.field28th = d(:,28);
    %tmp.field29th = d(:,29);

endfunction
