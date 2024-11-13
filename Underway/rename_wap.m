# script to rename DY151 WAP files that were in `non-standard' format

# In order to read the WAP files, the formats should be equivalent to: 
# 'DY151_2022061722_168_21_ASCII.013, where DY151 = cruise ID, 20220617 = YYYYMMDD, 22 = YY, 168 = jday,  21 = logger (2) and port (1) ID, ASCII = file/instrument type, .013 = WAP hour.

# This script standardises to above fromat for various special cases (uncomment each loop to run).

# wapdir = '/data/datasets/cruise_data/active/DY151/Data/Underway/WAP_Extracted/with_ACS167'
# wapdir = '/data/datasets/cruise_data/active/DY151/Data/Underway/WAP_Extracted/with_ACS167/'
wapdir = '/data/datasets/cruise_data/active/DY151/Data/Underway/WAP_Extracted/dark_counts_BB3/'    
    
# rename files with _v322_ extension
#files = glob(strcat(wapdir, 'DY151_*_v322_*'));
#for i=1:length(files)
#	oldname_i = strsplit(files{i}, '/'){10}
#	a = strsplit(oldname_i, '_'); # chunks of filename
#	if length(a) == 7
#		newname_i = strcat(a{1}, '_', a{2}, '22', '_', a{4}, '_', a{5}, '_', a{6}, '_', a{7})
#	elseif length(a) == 6
#		newname_i = strcat(a{1}, '_', a{2}, '22', '_', a{4}, '_', a{5}, '_', a{6})
#        elseif length(a) == 5
#		newname_i = strcat(a{1}, '_', a{2}, '22', '_', a{4}, '_', a{5})
#	endif
#	#rename(strcat(wapdir, oldname_i),strcat(wapdir, newname_i))
#endfor	

# rename files with _01.wlh22_ extension
#files = glob(strcat(wapdir, 'DY151_*_01.wlh22_*'));
#for i=1:length(files)
#	oldname_i = strsplit(files{i}, '/'){10}
#	a = strsplit(oldname_i, '_'); # chunks of filename
#	if length(a) == 7
#		newname_i = strcat(a{1}, '_', a{2}, '22', '_', a{4}, '_', a{5}, '_', a{6}, '_', a{7})
#	elseif length(a) == 6
#		newname_i = strcat(a{1}, '_', a{2}, '22', '_', a{4}, '_', a{5}, '_', a{6})
#       elseif length(a) == 5
#		newname_i = strcat(a{1}, '_', a{2}, '22', '_', a{4}, '_', a{5})#
#	endif
#	rename(strcat(wapdir, oldname_i),strcat(wapdir, newname_i))
#endfor	

# rename files with _v2.wlh22_ extension
#files = glob(strcat(wapdir, 'DY151_*_v2.wlh22_*'));
#for i=1:length(files)
#	oldname_i = strsplit(files{i}, '/'){10}
#	a = strsplit(oldname_i, '_'); # chunks of filename
#	if length(a) == 7
#		newname_i = strcat(a{1}, '_', a{2}, '22', '_', a{4}, '_', a{5}, '_', a{6}, '_', a{7})
#	elseif length(a) == 6
#		newname_i = strcat(a{1}, '_', a{2}, '22', '_', a{4}, '_', a{5}, '_', a{6})
#       elseif length(a) == 5
#		newname_i = strcat(a{1}, '_', a{2}, '22', '_', a{4}, '_', a{5})
#	endif
#	#rename(strcat(wapdir, oldname_i),strcat(wapdir, newname_i))
#endfor	

# rename files with _wlh22_ extension
#files = glob(strcat(wapdir, 'DY151_*.wlh22_*'));
#for i=1:length(files)
#	oldname_i = strsplit(files{i}, '/'){10}
#	a = strsplit(oldname_i, '_'); # chunks of filename
#	if length(a) == 7
#		newname_i = strcat(a{1}, '_', a{2}(1:8), '22', '_', a{4}, '_', a{5}, '_', a{6}, '_', a{7})
#	elseif length(a) == 6
#		newname_i = strcat(a{1}, '_', a{2}(1:8), '22', '_', a{4}, '_', a{5}, '_', a{6})
#       elseif length(a) == 5
#		newname_i = strcat(a{1}, '_', a{2}(1:8), '22', '_', a{4}, '_', a{5})
#	endif
#	rename(strcat(wapdir, oldname_i),strcat(wapdir, newname_i))
#endfor	


# rename files with _22_14222_ extension
#files = glob(strcat(wapdir, '*_22_14222_*'));
#for i=1:length(files)
#	oldname_i = strsplit(files{i}, '/'){10}
#	a = strsplit(oldname_i, '_'); # chunks of filename
#       newname_i = strcat(a{1}, '_', '2022052322_', a{3}, '_', a{4}, '_', a{5})
#	rename(strcat(wapdir, oldname_i),strcat(wapdir, newname_i))
#endfor	

# rename files with _22_14222_ extension
# files = glob(strcat(wapdir, '*_14222_*'));
# for i=1:length(files)
#	oldname_i = strsplit(files{i}, '/'){10}
#	a = strsplit(oldname_i, '_'); # chunks of filename
#       newname_i = strcat(a{1}, '_', a{2}, '_', a{4}, '_', a{5})
#	rename(strcat(wapdir, oldname_i),strcat(wapdir, newname_i))
#endfor	


# finally - rename strdate part of files using jday as ID (YYYYMMDD is often wrong!)
# inidate = "20220523"; # first get the correct mapping between jdays and strdates
# enddate = "20220623"; 
# [numdates, strdates, vecdates, jdays] = get_date_range(inidate, enddate);

#
# files = dir(strcat(wapdir,'*DY151*'));
# for i=1:length(files)
#	oldname_i = files(i).name;
#	a = strsplit(oldname_i, '_'); # chunks of filename
#	jday_i = str2num(a{3});
#	date_i = str2num(a{2});
#	for j = 1:length(jdays)
#		if jday_i == jdays(j)
#		   date_true = strcat(num2str(strdates(j,:)), '22');
#		   if strcmp(date_true, date_i) == 0
#		        oldname_i 
		        # jday_i
		        # date_i
		   	# date_true   
#			rest_of_id_i = [];
#			for k = 4: length(a)
#   				rest_of_id_i = strcat(rest_of_id_i, '_', a{k});
#			end
#			jday_i
#		   	newname_i = strcat(a{1}, '_', date_true, '_', num2str(jday_i), rest_of_id_i)
#		   	rename(strcat(wapdir, oldname_i),strcat(wapdir, newname_i))
#	           endif
#	       endif
#	endfor
#endfor

# with_ACS167
#inidate = "20220521"; # first get the correct mapping between jdays and strdates
#enddate = "20220523"; 
#[numdates, strdates, vecdates, jdays] = get_date_range(inidate, enddate);


#files = dir(strcat(wapdir,'*DY151*'));
#for i=1:length(files)
#	oldname_i = files(i).name
#	a = strsplit(oldname_i, '_'); # chunks of filename
#	jday_i = str2num(a{4});
#	for j = 1:length(jdays)
#		if jday_i == jdays(j)
#		   date_true = strcat(num2str(strdates(j,:)), '22');
#		   rest_of_id_i = [];
#		   for k = 5: length(a)
#			rest_of_id_i = strcat(rest_of_id_i, '_', a{k});
#		   end
#		   newname_i = strcat(a{1}, '_', date_true, '_', num2str(jday_i), rest_of_id_i)
#		   rename(strcat(wapdir, oldname_i),strcat(wapdir, newname_i))
#	       endif
#	    
#	endfor
#endfor


# dark_counts_BB3
#rename files with.wlh22_ extension
files = glob(strcat(wapdir, 'DY151_*.wlh22_*'))
#for i=1:length(files)
#	oldname_i = strsplit(files{i}, '/'){11}
	
#	a = strsplit(oldname_i, '_'); # chunks of filename
#	a{2} = a{2}(1:8);
#	if length(a) == 7
#		newname_i = strcat(a{1}, '_', a{2}, '22', '_',  a{3}, '_', a{4}, '_', a{5}, '_', a{6}, '_', a{7})
#	elseif length(a) == 6
#		newname_i = strcat(a{1}, '_', a{2}, '22', '_', a{3}, '_', a{4}, '_', a{5}, '_', a{6})
#       elseif length(a) == 5
#		newname_i = strcat(a{1}, '_', a{2}, '22', '_', a{3}, '_', a{4}, '_', a{5})
#	endif
#	rename(strcat(wapdir, oldname_i),strcat(wapdir, newname_i))
#endfor


#files = glob(strcat(wapdir, 'DY151_*_bb3_darks22_*'))
#for i=1:length(files)
#	oldname_i = strsplit(files{i}, '/'){11}
	
#	a = strsplit(oldname_i, '_'); # chunks of filename
#	if length(a) == 8
#		newname_i = strcat(a{1}, '_', a{2}, '22', '_',  a{5}, '_', a{6}, '_', a{7}, '_', a{8})
#	elseif length(a) == 7
#		newname_i = strcat(a{1}, '_', a{2}, '22', '_', a{5}, '_', a{6}, '_', a{7})
#       elseif length(a) == 6
#		newname_i = strcat(a{1}, '_', a{2}, '22', '_', a{5}, '_', a{6})
#	endif
#	rename(strcat(wapdir, oldname_i),strcat(wapdir, newname_i))
#endfor


