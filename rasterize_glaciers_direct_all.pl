#!/usr/bin/perl -w

#######################################################################
#
#	This code aggregates UAF PyGEM glacier model point data to raster files.
#	Suffix "direct" indicates no matching River network extents (no snapping).
#	This rasterized data is supposed to be used with high resolution river networks.
#
#	Written by Dr. A. Prusevich (alex.proussevitch@unh.edu)
#
#	August 2018
#		Last Modified-
#	Mar 2019	Added	Time series of glacier area.
#	Aug 2019	Added	Automation of init files.
#	Nov 2019	Added	Area scaling results/log file.
#	Feb 2020	Change	Time variable format. Now: YYYYMMDD as an integer.
#	Jan 2022	Change	Adopted to a new PyGEM output format. Other changes.
#
#######################################################################

use strict;
use File::Basename;
use File::Path;
use Geo::GDAL;
use Math::Trig qw/pi/;
use PDL;
use PDL::Char;
use PDL::Image2D;
use PDL::IO::FlexRaw;
use PDL::NetCDF;
use PDL::NiceSlice;
use Inline qw/Pdlpp/;
use Fcntl;
use Time::JulianDay;
use RIMS;
use RIMS::WBM;

my @time = (time());

use vars qw(*NEWERR *OLDERR);	# To avoid silly warning messages from GDAL, such as
open NEWERR, ">/dev/null";	# "No UNIDATA NC_GLOBAL:Conventions attribute"
open OLDERR, ">&STDERR";
STDOUT->autoflush(1);			# Disable buffering
$PDL::BIGPDL	= 1 if !$PDL::BIGPDL;	# Set BIGPDL, if needed (comment out if not needed)

set_autopthread_targ(4);	# Number of CPU threads
set_autopthread_size(1);	# Piddle size lower limit (in Meg) for multi-threading

#######################################################################
##################     Files and Directories        ###################

my $global_attrib = [
	['data_version',	'PyGEM data for PyGEM, Version 3 (January 2022)']];

my @model	= (
	'GFDL-ESM4',
	'CESM2',
	'MPI-ESM1-2-HR',
	'EC-Earth3',
	'NorESM2-MM',
	'INM-CM5-0',
	'MRI-ESM2-0',
	'BCC-CSM2-MR',
	'FGOALS-f3-L',
	'CESM2-WACCM',
	'EC-Earth3-Veg',
	'INM-CM4-8'
);
my @rcp		= qw(ssp126 ssp245 ssp370 ssp585);
my @RGI6_region	= map sprintf("%02d",$_), 1 .. 18;			# Use only regions needed for the Network grid
my $suffix	= 'MCMC_ba1_50sets_Global_glaciers';

my $data_dir	= '/net/nfs/saco/raid2/data/CMU_Glacier_Model_Output/globalsims_backup/simulations-cmip6/';
my $dir_out	= '/net/nfs/saco/raid2/data/CMU_Glacier_Model_Output/globalsims_backup/simulations-cmip6/rasterized/';
my $unzip_dir	= '/data/globalsims_backup/simulations-cmip6/';
my $netwrk_file	= '/net/nfs/swift/raid2/data/MERIT/MERIT_Hydro_IHU/15min/version_0.2/15min_flwdir_seg_endrh_v2.asc';

my $RGI6_area_file	= '/net/nfs/yukon/raid5/data/RGI6/global_rgi6_glaciers_1km.tif';

(my $grid_dir	= $dir_out . basename($netwrk_file)) =~ s/\.\w+$//;
 my @vars_m	= qw(runoff acc melt refreeze frontalablation massbaltotal prec);
 my @vars_y	= qw(area   volume);
 my $ref_year	= 2000;			# Reference year in the input data

			### Generic dataset metadata
my $meta = {'Var_Scale' => 1, 'Var_Offset' => 0, 'Processing' => '', 'Projection' => 'epsg:4326'};

#######################################################################
##################     Read Input Data        #########################

my $extent	= get_extent($netwrk_file);
my($lon, $lat)	= lonLat($extent);
my $grid	= [$lon, $lat, $$extent{gTransform}, $$extent{projection}];
my $cell_area	= cell_area($lon, $lat, $extent);					# Full cell area in km^2
my $RGI6_areaFr	= read_GDAL($extent,  $meta, 5, $RGI6_area_file, 1,0)->clip(0, 1);	# Resample by "average" (5)

#######################################################################
##################     Rasterize Glacier Data      ####################

printf "\nProcessing extends from %s\n", basename($netwrk_file);
(my $domain	= basename($netwrk_file)) =~ s/\.\w+$//;

foreach my $model (@model) {
foreach my $RCP (@rcp) {
  my $rcp =  $model eq 'ERA-Interim' ? '' : $RCP;
  print "\nWorking on Model (RCP): $model \($rcp\)\n";

	### Unzipping (modified script from Stanley Glidden
  unzip_data($data_dir, $unzip_dir, $model, $rcp, \@RGI6_region);

	### Make list of files
  my @file;
  foreach my $region (@RGI6_region) {
    my $dir_in	= "$unzip_dir/_unzipped/$model/$rcp/$region";
    opendir DIR,  $dir_in;
      push @file, map "$dir_in/$_", sort grep(m/\.nc$/, readdir(DIR));
    closedir DIR;
  }
  next unless @file;
 (my $file_out	= basename($file[0])) =~ s/.+($model.+)_all\.nc$/$grid_dir\/$model\/$rcp\/$1/;
  my $percent	= $#file/100;

	### Read dates
  my @t_time	= get_time_var($file[0], $ref_year);	# Julian Day
  my $t_ref	= julian_day(1900,1,1);
  my $t_month	= pdl(map($_ + 14		- $t_ref, @t_time));
  my $t_year	= pdl(map(julian_day($_, 7, 1)	- $t_ref, get_NetCDF_var($file[0],'year')->list));
  my @year	=					  get_NetCDF_var($file[0],'year')->list;

	### Prepare raster variables
  my %g_data_y	= map(($_ => zeroes($$extent{ncols}, $$extent{nrows}, $t_year ->dim(0))->copybad($$extent{mask})), @vars_y);
  my %g_data_m	= map(($_ => zeroes($$extent{ncols}, $$extent{nrows}, $t_month->dim(0))->copybad($$extent{mask})), @vars_m);
     $g_data_m{NN}	  =  zeroes($$extent{ncols}, $$extent{nrows})		       ->copybad($$extent{mask});

	### Read data
  for (my $i=0; $i <= $#file; $i++) {
    # printf("\r\tRasterizing %d of %d", $i+1, scalar(@file));
    printf("\r\tRasterizing %d %%", 100*$i/$#file) unless $i % $percent;
    my (%p_data_y, %p_data_m);

	### Read file data
    my $ncobj		= PDL::NetCDF->new ($file[$i], {MODE => O_RDONLY});
      my  $CenLon	= $ncobj->get('CenLon')->at(0);
      my  $CenLat	= $ncobj->get('CenLat')->at(0);
      map $p_data_m{$_}	= $ncobj->get('glac_'.$_.'_monthly'), @vars_m;
      map $p_data_y{$_}	= $ncobj->get('glac_'.$_.'_annual' ), @vars_y;
    $ncobj->close();

	### Rasterize
    my @colRow		= map int,Geo::GDAL::ApplyGeoTransform($$extent{igTransform},$CenLon,$CenLat);
		# Skip pixels outside the network domain
    next if $colRow[0] < 0 || $colRow[1] < 0 || $colRow[0] >= $$extent{ncols} || $colRow[1] >= $$extent{nrows};
		# Accumulate point data to the raster grids
    map $g_data_y{$_}->(@colRow,)->flat += $p_data_y{$_}->flat, @vars_y;
    map $g_data_m{$_}->(@colRow,)->flat += $p_data_m{$_}->flat, @vars_m;
	$g_data_m{NN}->(@colRow )       += 1;					# Glacier count in each pixel
  }
  $g_data_y{area}	*= 1e-6;				# Convert units to km2
  $g_data_y{volume}	*= 1e-9;				# Convert units to km3
  print "\r\tRasterizing 100 %                    \n";		# Use spaces to clear CR

	### Convert pixel accumulated glacial area km2 to spatially distributed area fractions
  my $g_area_mask;
  ($g_data_y{area_frac}, $g_area_mask) = area_frac($g_data_y{area}, $RGI6_areaFr, $cell_area, $$extent{mask}, \@year,
						"$grid_dir/$model\_$rcp.area_scl.log");

	### Save grids to NetCDF files
  print "\nSaving output files\n";
  my $file_out_y	= $file_out . '_y.nc';		unlink $file_out_y;
  my $file_out_m	= $file_out . '_m.nc';		unlink $file_out_m;
  my $file_out_c	= $file_out . '_c.tif';		unlink $file_out_c;
  my $file_out_a	= $file_out . '_a.tif';		unlink $file_out_a;
  my %dataToWrite_y	= (
    'area'	=> [$g_data_y{area},		'Glacier area',		'km2'],
    'area_frac'	=> [$g_data_y{area_frac},	'Glacier area fraction','fraction'],
    'volume'	=> [$g_data_y{volume},		'Glacier volume',	'km3']);
  my %dataToWrite_m	= (
    'runoff'	=> [$g_data_m{runoff},		'Glacier melt',	 	'm3/month'],
    'precip'	=> [$g_data_m{prec},		'Glacier precipitation','m3/month'],
    'snow'	=> [$g_data_m{acc},		'Glacier snowfall',	'm3/month'],
    'melt'	=> [$g_data_m{melt},		'Glacier meltwater',	'm3/month'],
    'refreezing'=> [$g_data_m{refreeze},	'Glacier refreezing',	'm3/month'],
    'ablation'	=> [$g_data_m{frontalablation},	'Frontal ablation',	'm3/month'],
    'massbaltot'=> [$g_data_m{massbaltotal},	'Mass balance total',	'm3/month']);
  write_nc( $file_out_y, 1, $t_year,  $grid, \%dataToWrite_y, {TS_RESOLUTION => 'yearly', ATTRIB => $global_attrib});
  write_nc( $file_out_m, 1, $t_month, $grid, \%dataToWrite_m, {TS_RESOLUTION => 'monthly',ATTRIB => $global_attrib});
  write_tif($extent,'Int16',$file_out_c, $g_data_m{NN});
  write_tif($extent,'Int16',$file_out_a, $g_area_mask);

	### Save init files
  my $iFile_runoff	= $grid_dir . "_init_files/$model\_$rcp\_$suffix\_glaciers_runoff_m.init";
  my $iFile_glMelt	= $grid_dir . "_init_files/$model\_$rcp\_$suffix\_glaciers_glIceMelt_m.init";
  my $iFile_volume	= $grid_dir . "_init_files/$model\_$rcp\_$suffix\_glaciers_volume_y.init";
  my $iFile_areaFr	= $grid_dir . "_init_files/$model\_$rcp\_$suffix\_glaciers_area_y.init";
  my $init_runoff	= init_runoff($file_out_m, $model, $rcp, $domain, $t_month, $t_ref);
  my $init_glMelt	= init_glMelt($file_out_m, $model, $rcp, $domain, $t_month, $t_ref);
  my $init_volume	= init_volume($file_out_y, $model, $rcp, $domain, $t_year,  $t_ref);
  my $init_areaFr	= init_areaFr($file_out_y, $model, $rcp, $domain, $t_year,  $t_ref);
  save_init($iFile_runoff, $init_runoff);
  save_init($iFile_glMelt, $init_glMelt);
  save_init($iFile_volume, $init_volume);
  save_init($iFile_areaFr, $init_areaFr);
  rmtree("$unzip_dir/_unzipped/$model");			# Remove the unzipped source files

  push @time, time();
  print  "$model\_$rcp is done!\n";
  printf "\nTime used - %d hours, %d minutes, and %d seconds\n", time_used($time[-2],$time[-1]);

  last unless $rcp;
}}

#######################################################################
						# Report Total Time
if ($#model + $#rcp) {
  printf "\nTotal time used - %d hours, %d minutes, and %d seconds\n",
	time_used($time[0],time());
  print "\n\nAll Done!\n\n";
}
close NEWERR;	close OLDERR;
exit;

#######################################################################
######################  Functions  ####################################

sub area_frac
{
  my ($data, $area, $cell_area, $network, $year, $log_file) = @_;
  print "\nCalculating area fraction\n";

	### Find location of glacier pixels
  my $gl_idx	= whichND($data->mv(2,0)->sumover > 0);
  unless ($gl_idx->dim(1)) {
    print "\tNo glaciers are within the Network.\n";
    return $data;
  }
	### Build masks of glacier catchments
  my $area_frac	= zeroes($data   ->dims)->copybad($network);
  my $mask	= zeroes($network->dims)->copybad($network);
  my $RGI6_rID	=($area > 0)->setbadtoval(0)->cc8compt->copybad($network); # Segmentation of RGI6 data to RGI mask IDs
  my $basins	= $RGI6_rID->max;
  foreach my $i ( 0 .. $gl_idx->dim(1)-1 ) {
    my @pixel	= $gl_idx(,$i)->list;
    unless ($mask->at(@pixel)) {
      my $id	= $RGI6_rID->at(@pixel);			# Use RGI6 segmentation IDs
         $id	= ++$basins	unless $id;			# Just in case DR point does not agree with RGI6
      my $IDX	= whichND($RGI6_rID == $id);
      my $idx	= whichND($network->upstreamMask(@pixel));	# Add upsteam area to the RGI6 mask
      $mask->indexND($idx) .= $id;
      $mask->indexND($IDX) .= $id;
  } }
  print "\tNumber of glacier catchments = $basins\n";

	## Scale glacier fractions within masks
  my @excessArea = ((0) x $data->dim(2));
  foreach my $basin (1 .. $basins) {
    print "\r\tProcessing catchment $basin";
# my $stop = 678;

    my $idx		= whichND($mask == $basin);
    my $basin_max	= 0;
    next unless $idx->dim(1);	# Multiple catchment IDs are overwritten by the highest ID in "mask" process above

    foreach my $i (0 .. $data->dim(2)-1) {
      my $trim		= 0;
      my $glc_area	=  $data(,,($i))->indexND($idx)->sum;
         $glc_area	=  $basin_max if $basin_max;
      my $prev_area	= ($i ? $area_frac(,,($i-1)) : $area)->indexND($idx) * $cell_area->indexND($idx);
      my $GLC_area	=  $prev_area->sum;
      unless ($GLC_area) {			# Case of zero previous glacial area: use RGI area to scale
         $prev_area	=  $area->indexND($idx) * $cell_area->indexND($idx);
         $GLC_area	=  $prev_area->sum;
      }
# printf "%s\t$i - Prev/Targ area = %f / %f\n", $i?'':"\n", $prev_area->sum, $glc_area if $basin == $stop;
		### Scaling !!!
      my($frac,$excess)	= scale($glc_area, $GLC_area, $prev_area, $cell_area->indexND($idx));
      $area_frac(,,($i))->indexND($idx) .= $frac;
# print "Excess $excess in catchment $basin \($i\)\n" if $excess && $basin == $stop;

		### Distribute excess area (if any) around border of glacier pixel within catchment
      while (abs($excess) > 1e-4) {
	my $b_mask	= ($mask == $basin);
# printf "b_mask = %d\n", $b_mask->sum if $basin == $stop;
	my $g_mask	=  $b_mask & ($area_frac(,,($i)) > 0);
	my $G_mask	=  $g_mask->setvaltobad(0)->patchbad2d->isgood & $b_mask;
	my $b_idx	=  whichND($G_mask - $g_mask);	# Index of border pixels

		### Overfilled catchment $basin. Expand catchment for the next downstream pixel
	unless ($b_idx->dim(1)) {
# print "Expanding!!!\n" if $basin == $stop;
	  my ($msk, $pixel, $fail) = grow_basin($b_mask, $mask, $network);
	  if ($fail) {	# Fail to expand
	    $basin_max		 = $cell_area->indexND($idx)->sum;
	    last;
	  } else {	# Good to expand
	    $trim		 = 1;
	    $idx		 = whichND($msk);
	    $mask->indexND($idx).= $basin;
	    $b_idx		 = pdl([$pixel]);
# 	      $b_idx = whichND($msk & (!$b_mask));
	} }
	$area_frac(,,($i))->indexND($b_idx) .= ($excess/$b_idx->dim(1) / $cell_area->indexND($b_idx))->hclip(1);
	$excess  -=  ($area_frac(,,($i))->indexND($b_idx) * $cell_area->indexND($b_idx))->sum;
# print  "\t\texcess = $excess : ",join('x',$b_idx->dims),"\n" if $basin == $stop;
      }
		### Trim overage catchment after expansion is finished
      if ($trim) {
	my $IDX			 = whichND(($mask == $basin) & ($area_frac(,,($i)) > 0));
	$mask->indexND($idx)	.= 0;
	$mask->indexND($IDX)	.= $basin;
	$idx			 = $IDX;
      }
      $excessArea[$i] += $excess;
    }
# die if $basin == $stop;
  }

	## Print scaling results
  prepare_dir($log_file);
  open (FILE,">$log_file") or die "Couldn't open $log_file, $!";
  print FILE	   "\t\t\t     Percent         Delta      Excess      Area\n";
  print "\n\tDone:\n\t\t\t     Percent         Delta      Excess      Area\n";
  foreach my $i (0 .. $data->dim(2)-1) {
    my $glc_area	= (	$area_frac(,,($i  ))		* $cell_area)->sum;
    my $prev_area	= (($i?	$area_frac(,,($i-1)) : $area)	* $cell_area)->sum;
    my $str		= sprintf "\t\tYear $$year[$i] = %7.2f %%  :%8.1f km2 :%6.1f km2 : %.1f km2\n",
			  $glc_area/$prev_area*100, $glc_area - $prev_area, $excessArea[$i], $glc_area;
    print FILE	$str;
    print	$str;
  }
  close FILE;

  return $area_frac, $mask;
}

sub scale
{
  my ($t_sum, $p_sum, $area, $cell_area) = @_;
  my  $excess	= 0;
		### Decrease area
  if ($t_sum <= $p_sum) {
    $area *= $t_sum / $p_sum;
  }
		### Increase area
  else {
    while(abs($p_sum - $t_sum)/$t_sum > 1e-4) {
      my $idx	= whichND(($area > 0) & ($area < $cell_area));
      unless ($idx->nelem) { $excess = $t_sum - $p_sum; last; }

      my $sum	= $area->indexND($idx)->sum;
      my $scale	=($sum + ($t_sum - $p_sum)) / $sum;
      $area->indexND($idx) .= ($area->indexND($idx) * $scale)->hclip($cell_area->indexND($idx));
      $p_sum 	= $area->sum;
    }
  }
  return $area / $cell_area, $excess;
}

sub grow_basin
{
  my($mask, $MASK, $network)	= @_;
  my @point		= whichND($mask)->(,0)->list;
  my $basin		= $MASK->at(@point);
  my @direction		= ([0,1,1,0,-1,-1,-1,0,1],[0,0,1,1,1,0,-1,-1,-1]);
  my %log2		= (0=>0,1=>1,2=>2,4=>3,8=>4,16=>5,32=>6,64=>7,128=>8);
  my $fail		=  0;

  while ($MASK->at(@point)) {		# It is possibel to step over previously expanded mask
    my $dir	= $log2{$network->at(@point)};
    my @pnt	= map	$point[$_]+$direction[$_][$dir], 0..1;
       $fail	= 1 if	$pnt[0]<0 || $pnt[1]<0 || $pnt[0]==$network->dim(0) || $pnt[1]==$network->dim(1) ||
			$network(@pnt)->isbad  || $dir==0;
       @point	= @pnt	unless	$fail;
    last		if	$fail;
  }
  my $add	= $network->upstreamMask(@point) - $mask;
     $add->where($MASK) .= 0;		# Do not use expanded upsteam area with other catchment IDs

  return $mask+$add, \@point, $fail;
}
#######################################################################

sub unzip_data
{
  my ($base, $base2, $GCM, $ssp, $regions) = @_;
  my  $inpathfile = $base . join('_', $GCM, $ssp) . '_stats.zip';

  foreach my $region (@$regions) {
    printf "\r\tUnzipping region $region of %d", scalar(@$regions);
    my $inpath     = $base   .  '_zipped/' . "$region/stats/";
    my $outpath    = $base2  .'_unzipped/' . "$GCM/$ssp/$region/";
    my $inpathfile = $inpath . join('_', $GCM, $ssp) . '_stats.zip';
	# Check/make output directory
    next if -d $outpath;
    mkpath($outpath,0,0775) or die "Cannot make directory...\n$outpath\n";
	# Unzip and change permissions
    system "unzip -qq -d $outpath $inpathfile";
    system "find $outpath -type d | xargs chmod 775";
    system "find $outpath -type f | xargs chmod 664";
  } print "\n";
}

#######################################################################

sub get_NetCDF_var
{
  my ($file,$var_name) = @_;

  my $ncobj = PDL::NetCDF->new ($file, {MODE => O_RDONLY});
    my $var = $ncobj->get($var_name);
  $ncobj->close();

  return $var;
}

#######################################################################

sub get_time_var
{
  my($file, $ref_year)	= @_;
  my @data	= get_NetCDF_var($file,'time')->list;
  my $t_ref	= julian_day($ref_year,1,1);
  my @j_date	= map $_ + $t_ref, @data;

  return @j_date;
}

#######################################################################

sub get_NetCDF_textvar
{
  my ($file,$var_name) = @_;

  my  @text;
  my  $text = `ncks -H -v $var_name $file`;
	# Different versions of "ncks"
  if ($text =~ m/^netcdf/) {			### Version on magma
   my $str  = $text=~m/data:\s+glac_attrs = (.+) ;/ ? $1 : die "Error in \"get_NetCDF_textvar\" function...\n";
      $str  =~ s/"//g;
      @text = split ', ', $str;
  }
  else {					### Version on yukon
      @text = split m/\s*\n/, $text;
      map s/.+?=//, @text;
  }

  return @text;
}

#######################################################################

sub get_NetCDF_attr
{
  my ($file,$var_name,$att_name) = @_;

  my $ncobj = PDL::NetCDF->new ($file, {MODE => O_RDONLY});
    my $att = $ncobj->getatt($att_name,$var_name);
  $ncobj->close();

  return $att;
}

#######################################################################

sub init_runoff
{
  my ($file, $model, $rcp, $domain, $time, $t_ref) = @_;
  my  @date	= $time->list;
  my  $dateBeg	= sprintf "%04d-%02d-00", (inverse_julian_day($date[ 0] + $t_ref))[0,1];
  my  $dateEnd	= sprintf "%04d-%02d-00", (inverse_julian_day($date[-1] + $t_ref))[0,1];
  my  $bands	= scalar   @date;

  my $str = <<EOF;
{
  Code_Name	=> '$model\_$rcp\_$suffix\_glaciers_runoff_m',
  Time_Series	=> 'monthly',
  Start_Date	=> '$dateBeg',
  End_Date	=> '$dateEnd',
  Name		=> 'PyGEM $model $rcp RGI-5 glacier Runoff rasterized to $domain grid (UAF-Rounce v.3), monthly',
  Var_Name	=> 'runoff',
  Units		=> 'm3/sec',
  Var_Scale	=> 3.858e-07,
  Var_Offset	=> 0,
  Orig_Units	=> 'm3/month',
  Bands		=> $bands,
  Processing	=> 'noScale',
  Projection	=> 'epsg:4326',
  File_Path	=> '(1e20,1e5):$file;'
}
EOF

  return $str;
}

sub init_glMelt
{
  my ($file, $model, $rcp, $domain, $time, $t_ref) = @_;
  my  @date	= $time->list;
  my  $dateBeg	= sprintf "%04d-%02d-00", (inverse_julian_day($date[ 0] + $t_ref))[0,1];
  my  $dateEnd	= sprintf "%04d-%02d-00", (inverse_julian_day($date[-1] + $t_ref))[0,1];
  my  $bands	= scalar   @date;

  my $str = <<EOF;
{
  Code_Name	=> '$model\_$rcp\_$suffix\_glaciers_glIceMelt_m',
  Time_Series	=> 'monthly',
  Start_Date	=> '$dateBeg',
  End_Date	=> '$dateEnd',
  Name		=> 'PyGEM $model $rcp RGI-5 glacier Runoff rasterized to $domain grid (UAF-Rounce v.3), monthly',
  Var_Name	=> 'melt',
  Units		=> 'm3/sec',
  Var_Scale	=> 3.858e-07,
  Var_Offset	=> 0,
  Orig_Units	=> 'm3/month',
  Bands		=> $bands,
  Processing	=> 'noScale',
  Projection	=> 'epsg:4326',
  File_Path	=> '(1e20,1e5):$file;'
}
EOF

  return $str;
}

sub init_volume
{
  my ($file, $model, $rcp, $domain, $time, $t_ref) = @_;
  my  @date	= $time->list;
  my  $dateBeg	= sprintf "%04d-00-00", (inverse_julian_day($date[ 0] + $t_ref))[0];
  my  $dateEnd	= sprintf "%04d-00-00", (inverse_julian_day($date[-1] + $t_ref))[0];
  my  $bands	= scalar   @date;

  my $str = <<EOF;
{
  Code_Name	=> '$model\_$rcp\_$suffix\_glaciers_volume_y',
  Time_Series	=> 'yearly',
  Start_Date	=> '$dateBeg',
  End_Date	=> '$dateEnd',
  Name		=> 'PyGEM $model $rcp RGI-5 glacier Volume rasterized to $domain grid (UAF-Rounce v.3), yearly',
  Var_Name	=> 'volume',
  Units		=> 'm3',
  Var_Scale	=> 1e9,
  Var_Offset	=> 0,
  Orig_Units	=> 'km3',
  Bands		=> $bands,
  Processing	=> 'noScale',
  Projection	=> 'epsg:4326',
  File_Path	=> '(1e20,1e5):$file;'
}
EOF

  return $str;
}

sub init_areaFr
{
  my ($file, $model, $rcp, $domain, $time, $t_ref) = @_;
  my  @date	= $time->list;
  my  $dateBeg	= sprintf "%04d-00-00", (inverse_julian_day($date[ 0] + $t_ref))[0];
  my  $dateEnd	= sprintf "%04d-00-00", (inverse_julian_day($date[-1] + $t_ref))[0];
  my  $bands	= scalar   @date;

  my $str = <<EOF;
{
  Code_Name	=> '$model\_$rcp\_$suffix\_glaciers_area_y',
  Time_Series	=> 'yearly',
  Start_Date	=> '$dateBeg',
  End_Date	=> '$dateEnd',
  Name		=> 'PyGEM $model $rcp RGI-5 glacier area fraction rasterized to $domain grid (UAF-Rounce v.3), yearly',
  Var_Name	=> 'area_frac',
  Units		=> 'fraction',
  Var_Scale	=> 1,
  Var_Offset	=> 0,
  Orig_Units	=> 'fraction',
  Bands		=> $bands,
  Processing	=> 'noScale',
  Projection	=> 'epsg:4326',
  File_Path	=> '(1e20,1e5):$file;'
}
EOF

  return $str;
}

#######################################################################

sub save_init
{
  my ($file,$str) = @_;
  prepare_dir($file);
  open (FILE,">$file") or die "Couldn't open $file, $!";
    print FILE $str;
  close FILE;
}

#######################################################################
###################  PDL::PP Functions  ###############################

__DATA__
__Pdlpp__

#######################################################################

pp_addhdr('
  #include <unistd.h>       /* we need defs of XXXX */
  #include <stdio.h>

  static void upstr(PDL_Byte * flowDir, PDL_Long * stack, PDL_Long n_size, PDL_Long m_size, PDL_Long N, PDL_Long M)
  {
    long from[2][8] = { {1,1,0,-1,-1,-1,0,1} , {0,1,1,1,0,-1,-1,-1} };
    long dir, ind, xx, yy;
    long NN = N;	long count = 1;
    long MM = M;	long pos   = 0;
    stack[0] = N + M*n_size;

    while (pos < count) {
      MM = stack[pos] / n_size;
      NN = stack[pos] - n_size*MM;
      pos++;
      for (dir=0; dir<8; dir++) {
        xx  = NN - from[0][dir];
        yy  = MM - from[1][dir];
        if (xx<0 || yy<0 || xx==n_size || yy==m_size) break;	// This never should happen
        ind = xx + yy*n_size;
        if (flowDir[ind] == (0x01<<dir)) stack[count++] = ind;
      }
    }
  }
');

#######################################################################

pp_def('upstreamMask', HandleBad => 1,
  Pars => 'byte flowDir(n,m);
    int N(); int M();
    byte [o] mask(n,m);',
  Code => '
    int ind;
    int n_size = $SIZE(n);    int NN = $N();
    int m_size = $SIZE(m);    int MM = $M();
    int *myStack;	myStack = malloc(n_size*m_size*sizeof *myStack);

    loop(n,m) %{ $mask() = 0; %}  	//	Initialization of the output arrays
    for (ind=0; ind < n_size*m_size; ind++) {
      myStack[ind] = -1;
    }

    upstr($P(flowDir),myStack,n_size,m_size,NN,MM);
    ind = 0;
    while (myStack[ind] != -1 && ind < n_size*m_size) {
      MM = myStack[ind] / n_size;
      NN = myStack[ind] - n_size*MM;
      ind++;
      $mask(n=>NN,m=>MM) = 1;
    }
	// Free OS memory
    free(myStack);
');

#######################################################################

pp_def('upstrAccumAll', HandleBad => 1,
  Pars => 'byte flowDir(n,m); double data(n,m);
    double [o] upstrAccData(n,m);',
  Code => '
    int i, j, ind;
    int n_size = $SIZE(n);
    int m_size = $SIZE(m);
    double *myData;	myData  = malloc(n_size*m_size*sizeof *myData);
    int    *myStack;	myStack = malloc(n_size*m_size*sizeof *myStack);

	//	Initialization of the output arrays
    loop(n,m) %{ $upstrAccData() = 0; %}
    for (j=0; j<m_size; j++) {
      for (i=0; i<n_size; i++) {
	ind = i + j*n_size;
	myData[ind]  = ( $ISBAD($data(n=>i,m=>j)) ) ? 0 : $data(n=>i,m=>j);
	myStack[ind] = -1;
      }
    }

    for (j=0; j<m_size; j++) {
      for (i=0; i<n_size; i++) {
	if ( $ISBAD($flowDir(n=>i,m=>j)) ) {
	  $SETBAD($upstrAccData(n=>i,m=>j));
	}
	else {
	  upstr($P(flowDir),myStack,n_size,m_size,i,j);
	  ind = 0;
	  while (myStack[ind] != -1 && ind < n_size*m_size) {
	    $upstrAccData(n=>i,m=>j) += myData[myStack[ind]];
	    myStack[ind++] = -1;				// Reset myStack
	  }
	}
      }
    }
	// Free OS memory
    free(myData);
    free(myStack);
');

#######################################################################

pp_done();
