%function tmp = step2_uway();
   % Compute bio-physical quantities (IOPs) from (uncallibrated) optical data

   clear all
   PLOT = 1;

   % Load paths and common variables
   #run("input_parameters.m")
   run("../input_parameters.m")


   global OUT_PROC # this is for plot_spectra2.m
   global UWAY_DIR

   % Create date range
   [numdates, strdates, vecdates, jday_in] = get_date_range(inidate,enddate);

   fn_saved = glob([DIR_STEP1 "*mat"]);

   global YYYY = vecdates(1, 1); % Assumes all AMT days are within same year!! % used as processing Id

   % Change first day to process more than just last sampled day
   first_day = find_index_strdate_in_glob(fn_saved, sprintf("%d", jday_in(1))); % follows from ini and end dates
   last_day = find_index_strdate_in_glob(fn_saved, sprintf("%d", jday_in(end)));

   % Need to overwrite array of jdays with dates from saved files
   for ifile = 1:length(fn_saved)
       jdays(ifile) = str2num(strsplit(fn_saved{ifile}, "."){1}(end-2:end)); % creates jday array-define
   endfor

   dailyfiles = dir([DIR_STEP1 "*mat"]  ); % redundancy with line 27? (just different format)


   for iday = last_day -2 : last_day

        disp(["\n---------" dailyfiles(iday).name "--------\n"] )
        fflush(stdout);


       # initialize output structure with nans
           ini_out(dailyfiles(iday).name, jdays(iday));

        % First process Ship ctd data
        % (needed by bb3 processing)
        disp("\nprocessing SHIPs UNDERWAY data...");



	     % Cook version of underway function (AMT31)
        if jdays(iday) > 354
        uway = step2h_underway_cook_make_processed(jdays(iday), strdates(iday,:),
                FUNC_GGA,
                DIR_GPS, FN_GPS,
                DIR_ATT, FN_ATT,
                DIR_DEPTH, FN_DEPTH,
                DIR_TS, FN_SURF, FN_METDATA, FN_LIGHT,
                DIR_TSG, FN_TSG);
        elseif jdays(iday) < 354
                uway = step2h_underway_cook_make_processed(jdays(iday), strdates(iday,:),
                FUNC_GGA,
                DIR_GPS, FN_GPS,
                DIR_ATT, FN_ATT,
                DIR_DEPTH, FN_DEPTH,
                DIR_TS, FN_SURF, FN_METDATA, FN_LIGHT,
                DIR_TSG, FN_TSG);
        elseif  jdays(iday) == 354
                uway = step2h_underway_cook_make_processed(jdays(iday), strdates(iday,:),
                FUNC_GGA,
                DIR_GPS, FN_GPS,
                DIR_ATT, FN_ATT,
                DIR_DEPTH, FN_DEPTH,
                DIR_TS, FN_SURF, FN_METDATA, FN_LIGHT,
                DIR_TSG, FN_TSG);
        disp("day with multiple files - need to check output");
         end

        disp("...done");

        jday_str = dailyfiles(iday).name(end-6:end-4);

        % Load WAPvars from step1 output file
        load([DIR_STEP1 dailyfiles(iday).name]);


        % Idea is that flow is always there
        % (also needed by ac9 processing)
        disp("processing Flow data...");
        flow = step2f_flow_make_processed(WAPvars.flow, dailyfiles(iday));
        disp("...done");

        % Cycle through the variables within WAPvars
        instruments = fieldnames(WAPvars);
        for iWAP = 1:length(instruments)

           disp(["Processing ", instruments{iWAP}, " data..."]);

           switch instruments{iWAP}
               case "flow"
                   disp("Flow already processed")

               case "acs"
                   if  sum(~isnan(WAPvars.acs.raw.med(:,1))) > 1;
                       step2a_acs_amt_make_processed(WAPvars.acs, dailyfiles(iday), iday, acs_lim, FORCE=0, "acs");
                   endif

               case "acs2"
                   if  sum(~isnan(WAPvars.acs2.raw.med(:,1)))> 1
                       step2a_acs_amt_make_processed(WAPvars.acs2, dailyfiles(iday), iday, acs_lim, FORCE=0, "acs2");
                   endif

               case "ac9"
                   step2a_ac9_amt_make_processed(WAPvars.ac9, dailyfiles(iday), ac9_lim, FORCE=0, flow);

## uncomment this when you want to process BB3 data
#               case "bb3"
#                   step2b_bb3_amt_make_processed(WAPvars.bb3, uway, dailyfiles(iday), iday, bb_opt_lim, CRUISE);

               case "cstar"
                   step2d_cstar_make_processed(WAPvars.cstar, dailyfiles(iday), cstar_lim);

               case "ctd"
                   step2f_ctd_make_processed(WAPvars.ctd, dailyfiles(iday));

#               otherwise
#                   disp("Instrument to be implemented")
#                   keyboard

           endswitch
           disp("...done");
       endfor
       disp("\n");
       toc
   endfor

 #  keyboard


 if PLOT == 1
   % Plot spectra from acs
   disp("\nplotting spectra...");
   for iday = first_day:last_day
#      disp(num2str(jdays(iday)));
#      fflush(stdout);
      plot_spectra2(sprintf("%d",jdays(iday)),spectra_alim, spectra_clim, chl_lim);
   endfor
 endif

   %% save chl for Bob
   %    secchi_chl
%endfunction




