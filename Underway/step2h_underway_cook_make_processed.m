function tmp = step2h_underway_cook_make_processed(doy, date_str, FUNC_GGA, DIR_GPS, FN_GPS, DIR_ATT, FN_ATT, DIR_DEPTH, FN_DEPTH, DIR_TS, FN_SURF, FN_METDATA, FN_LIGHT, DIR_TSG, FN_TSG)


   # This function was developed by tjor in August 2022, so that discovery underway metadata could be read/processed in a simlar way to
   # JCR data for back-processing.

   # The function combines `step2h_underway_make_processed.m' (2022 version of underway processing file,
   # which acts on a single day at a time) with `step2h_underway_amt27_makeprocessed.m'(version which contains
   # file reading syntax & functions for the discovery)

   # INPUTS:  doy = julian day and date_str = YYYYMMDD format. All other inputs are specified in inputparamters.m, and
   # remove previous hardcoding within the discovery rd functions.

   % Global variables from step2
   global YYYY
   global DIR_STEP1
   global FN_ROOT_STEP2


   % Filenames for meta data
   fn_gps = glob([DIR_GPS date_str FN_GPS]); % used in gga function
   fn_att = glob([DIR_ATT date_str FN_ATT]);
   fn_depth = glob([DIR_DEPTH date_str FN_DEPTH]); %

   fn_surf = glob([DIR_TS date_str FN_SURF]); % used in oceanlogger function
   fn_met = glob([DIR_TS date_str FN_METDATA]); %
   fn_light = glob([DIR_TS date_str FN_LIGHT]); %
   fn_tsg = glob([DIR_TSG date_str FN_TSG]); %

   disp(fn_gps{1})
   disp(fn_depth{1})
   disp(fn_att{1})

   disp(fn_surf{1})
   disp(fn_met{1})
   disp(fn_light{1})
   disp(fn_tsg{1})

 #  keyboard
   disp('Processing ship''s underway data...')


   # 20/12 - temp
  # if doy ==355
  #  fn_gps = fn_gps(2);
  #  fn_depth = fn_depth(2);
  #  fn_att = fn_att(2);

   # fn_surf = fn_surf(2);
   # fn_met = fn_met(2);
   # fn_light = fn_light(2);
    #fn_tsg = fn_tsg(2);
 #  end

  # if doy ==354
  #     fn_gps = fn_gps(1:2);
  #     fn_depth = fn_depth(1:2);
  #     fn_att = fn_att(1:2);

  #     fn_surf = fn_surf(1:2);
  #     fn_met = fn_met(1:2);
  #     fn_light = fn_light(1:2);
  #     fn_tsg = fn_tsg(1:2);
 #  endif

   % load meta data and gps files - files are passed as arguments
   tmp1 = rd_seatech_gga_cook(fn_gps, fn_att, fn_depth); # file reader same on disco and cook
   tmp2 = rd_oceanlogger_cook(fn_surf, fn_met, fn_light, fn_tsg); # cook sometimes has multiple files


   % create daily time vector with one record per minute of the day (24*60=1440) - note: lines 40-
   tmp.time = y0(YYYY)-1 + doy + [0:1440-1]'/1440; # time vector to match 1-min binned optics data

  # keyboard

    %interpolate underway data to one-minute samples (and combine in single data structure)
    flds1 = fieldnames(tmp1);
    for ifld1=2:length(flds1) % index 2 - skips time field
         tmp.(flds1{ifld1}) = nan(size(tmp.time));
         if ~isempty(tmp1.time)
            tmp.(flds1{ifld1}) = interp1(tmp1.time, tmp1.(flds1{ifld1}), tmp.time,'extrap')
         endif
    endfor


    flds2 = fieldnames(tmp2);
    for ifld2=2:length(flds2) % skips time field
         tmp.(flds2{ifld2}) = nan(size(tmp.time));
         if ~isempty(tmp2.time)
                  flds2(ifld2)
            tmp.(flds2{ifld2}) = interp1(tmp2.time, tmp2.(flds2{ifld2}), tmp.time,'extrap');
         endif
    endfor

    % save underway ship's data to optics file
    savefile = [FN_ROOT_STEP2 num2str(doy) '.mat'];
    if (exist(savefile))
         load(savefile);
    endif

    out.uway = tmp;
    save('-v6', savefile , 'out' );


endfunction
