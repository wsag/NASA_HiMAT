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
#
#######################################################################

use strict;
use File::Basename;
use File::Path;
use Geo::GDAL;
use Geo::Proj4;
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

my ($init_file, $io_file) = map '/net/home/cv/alexp/perl/wbm/wbm_dev/'.$_, qw(wbm_path.init wbm_io.pl);
{			### wbm_io.pl must be in the same directory
  local @ARGV = ($init_file);
  require $io_file;
}
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
	['data_version',	'PyGEM data for HiMAT, Version 3 (September 2019)']];

my @model	= qw(	CCSM4);		my @rcp = (85);
# my @model	= qw(	ERA-Interim);	my @rcp = (85);
# my @model	= qw(	NorESM1-M	MPI-ESM-LR	MPI-ESM-MR	HadGEM2-ES	FGOALS-g2
# 			CNRM-CM5	CanESM2
# 			NorESM1-ME	MRI-CGCM3	MIROC-ESM	MIROC-ESM-CHEM	MIROC5
# 			IPSL-CM5A-LR	IPSL-CM5A-MR	GFDL-CM3	GFDL-ESM2M	GFDL-ESM2G
# 			GISS-E2-R	CSIRO-Mk3-6-0	CESM1-CAM5	CCSM4		bcc-csm1-1);
# my @rcp		= (26, 45, 60, 85);

my $data_dir	= '/net/nfs/merrimack/raid2/data/glaciers_6.0/';
# my $netwrk_file	= '/net/nfs/zero/home/WBM_TrANS/data/nepal_1km_v2.asc';
# my $netwrk_file	= '/net/nfs/zero/home/WBM_TrANS/data/karsub3.asc';
my $netwrk_file	= '/net/nfs/zero/home/WBM_TrANS/data/HiMAT_full_210_Subset.asc';

my $RGI6_area_file	= '/net/nfs/yukon/raid5/data/RGI6/global_rgi6_glaciers_1km.tif';

(my $grid_dir	= $data_dir . basename($netwrk_file)) =~ s/\.\w+$//;
 my @vars_m	= qw(runoff acc melt refreeze frontalablation massbaltotal prec);
 my @vars_y	= qw(area   volume);

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
  my $rcp =  $model eq 'ERA-Interim' ? '' : 'rcp'.$RCP;
  print "\nWorking on Model (RCP): $model \($rcp\)\n";

	### Make list of files
#   opendir DIR,	      $data_dir. 'trishuli_naltar_v3' or die "Cannot open directory $data_dir\trishuli_naltar_v3/";
#     my @file 	= map $data_dir. "trishuli_naltar_v3/$_", sort grep(m/.+$model\_$rcp.+2017\.nc$/, readdir(DIR));
  opendir DIR,	      $data_dir. $model;
    my @file 	= map $data_dir."$model/$_", sort grep(m/.+$rcp.+\.nc$/, readdir(DIR));
  closedir DIR;

  next unless @file;
#   ( my $file_out= basename($file[0])) =~ s/^R\d{6}_(.+)\.nc$/$grid_dir\/$1/;
#   ( my $file_out= basename($file[0])) =~ s/^R\d{2}_(.+)\.nc$/$grid_dir\/$1/;
  ( my $file_out= basename($file[0])) =~ s/^R\d{2}--all--(.+)\.nc$/$grid_dir\/$1/;

# 	### Read dates in format:	time:units = "days since 1999-10-01"
#   my $t_att_m	= get_NetCDF_attr($file[0],'time','units');	$t_att_m = $1 if $t_att_m =~ m/(\d{4}-\d{2}-\d{2})/;
#   die "Failed to read reference time. Aborting...\n" unless			 $t_att_m =~ m/^\d{4}-\d{2}-\d{2}$/;
#   my $t_ref_m	= julian_day(split m/-/, $t_att_m);
#   my $t_ref	= julian_day(1900,1,1);
#   my $t_year	= pdl(map(julian_day($_, 7, 1)	- $t_ref, get_NetCDF_var($file[0],'year_plus1')	->list));
#   my $t_month	= pdl(map($_ + 14 + $t_ref_m	- $t_ref, get_NetCDF_var($file[0],'time')	->list));
#   my @year	=					  get_NetCDF_var($file[0],'year_plus1')	->list;

	### Read dates in format:	YYYYMMDD as an integer
  my @t_time	= get_time_var($file[0]);	# Forcing the reference date in it to 1900-01-01
  my $t_ref_m	= julian_day(1900,1,1);
  my $t_ref	= julian_day(1900,1,1);
  my $t_month	= pdl(map($_ + 14 + $t_ref_m	- $t_ref, @t_time));
  my $t_year	= pdl(map(julian_day($_, 7, 1)	- $t_ref, get_NetCDF_var($file[0],'year_plus1')	->list));
  my @year	=					  get_NetCDF_var($file[0],'year_plus1')	->list;

	### Prepare raster variables
  my %g_data_y	= map(($_ => zeroes($$extent{ncols}, $$extent{nrows}, $t_year ->dim(0))->copybad($$extent{mask})), @vars_y);
  my %g_data_m	= map(($_ => zeroes($$extent{ncols}, $$extent{nrows}, $t_month->dim(0))->copybad($$extent{mask})), @vars_m);
     $g_data_m{NN}	  =  zeroes($$extent{ncols}, $$extent{nrows})		       ->copybad($$extent{mask});

	### Read data
  foreach my $data_file (@file) {
    printf "\nProcessing file- %s\n", $data_file;
    print  "\tReading data\n";
    my (%p_data_y, %p_data_m);

	### Read metadata
    my $gl_table	= get_NetCDF_var(    $data_file,'glacier_table');
    my @header		= get_NetCDF_textvar($data_file,'glac_attrs');
    my %hdr		= map(($header[$_] => $_), 0 .. $#header);
    my $coordCol	= pdl($hdr{CenLon}, $hdr{CenLat});	# X for (Lon, Lat) in the $gl_table
    my $percent		= $gl_table->dim(1)/100;

	### Read annual point data
    map $p_data_y{$_} = get_NetCDF_var($data_file, $_.'_glac_annual' )->((0),,), @vars_y;

	### Read monthly point data
    foreach my $var (@vars_m) {
	  $p_data_m{$var}         = get_NetCDF_var($data_file, $var.'_glac_monthly')->((0),,);	next if $var eq 'runoff';
	### Convert units from m to m3/month and use hydrological year shift as 3 below
      map $p_data_m{$var}->($_,) *= $p_data_y{area}->(int(($_+3)/12),) * 1e6, 0 .. $p_data_m{$var}->dim(0)-1;
    }
	### Rasterize
    for (my $i=0; $i < $gl_table->dim(1); $i++) {
      printf("\r\tRasterizing %d %%", 100*$i/$gl_table->dim(1)) unless $i % $percent;

      my @colRow	= map int,Geo::GDAL::ApplyGeoTransform($$extent{igTransform},$gl_table($coordCol,$i)->list);
		# Skip pixels outside the network domain
      next if $colRow[0] < 0 || $colRow[1] < 0 || $colRow[0] >= $$extent{ncols} || $colRow[1] >= $$extent{nrows};
		# Accumulate point data to the raster grids
      map $g_data_y{$_}->(@colRow,)->flat += $p_data_y{$_}->(,$i)->flat, @vars_y;
      map $g_data_m{$_}->(@colRow,)->flat += $p_data_m{$_}->(,$i)->flat, @vars_m;
	  $g_data_m{NN}->(@colRow )       += 1;					# Glacier count in each pixel
    }
    print "\r\tRasterizing 100 %\n";
  }
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
  my $iFile_runoff	= $grid_dir . "_init_files/$model\_$rcp\_c2_HiMAT_glaciers_runoff_m.init";
  my $iFile_glMelt	= $grid_dir . "_init_files/$model\_$rcp\_c2_HiMAT_glaciers_glIceMelt_m.init";
  my $iFile_volume	= $grid_dir . "_init_files/$model\_$rcp\_c2_HiMAT_glaciers_volume_y.init";
  my $iFile_areaFr	= $grid_dir . "_init_files/$model\_$rcp\_c2_HiMAT_glaciers_area_y.init";
  my $init_runoff	= init_runoff($file_out_m, $model, $rcp, $domain, $t_month, $t_ref);
  my $init_glMelt	= init_glMelt($file_out_m, $model, $rcp, $domain, $t_month, $t_ref);
  my $init_volume	= init_volume($file_out_y, $model, $rcp, $domain, $t_year,  $t_ref);
  my $init_areaFr	= init_areaFr($file_out_y, $model, $rcp, $domain, $t_year,  $t_ref);
  save_init($iFile_runoff, $init_runoff);
  save_init($iFile_glMelt, $init_glMelt);
  save_init($iFile_volume, $init_volume);
  save_init($iFile_areaFr, $init_areaFr);

  push @time, time();
  print  "$rcp is done!\n";
  printf "\nTime used - %d hours, %d minutes, and %d seconds\n", time_used($time[-2],$time[-1]);

  last unless $rcp;
}}

#######################################################################
						# Report Total Time
printf "\nTime used - %d hours, %d minutes, and %d seconds\n",
	time_used($time[0],time());
print "\n\nAll Done!\n\n";

close NEWERR;	close OLDERR;
exit;

#######################################################################
######################  Functions  ###################################get_NetCDF_var#

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
  my $file	= shift;
  my @data	= get_NetCDF_var($file,'time')->list;
  my $t_ref	= julian_day(1900,1,1);
  my @days	= map julian_day(substr($_,0,4), substr($_,4,2), substr($_,6,2)) - $t_ref, @data;

  return @days;
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
  Code_Name	=> '$model\_$rcp\_c2_HiMAT_glaciers_runoff_m',
  Time_Series	=> 'monthly',
  Start_Date	=> '$dateBeg',
  End_Date	=> '$dateEnd',
  Name		=> 'HiMAT $model $rcp RGI-6.0 glacier Runoff rasterized to $domain grid (UAF-Rounce v.3), monthly',
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
  Code_Name	=> '$model\_$rcp\_c2_HiMAT_glaciers_glIceMelt_m',
  Time_Series	=> 'monthly',
  Start_Date	=> '$dateBeg',
  End_Date	=> '$dateEnd',
  Name		=> 'HiMAT $model $rcp RGI-6.0 glacier Runoff rasterized to $domain grid (UAF-Rounce v.3), monthly',
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
  Code_Name	=> '$model\_$rcp\_c2_HiMAT_glaciers_volume_y',
  Time_Series	=> 'yearly',
  Start_Date	=> '$dateBeg',
  End_Date	=> '$dateEnd',
  Name		=> 'HiMAT $model $rcp RGI-6.0 glacier Volume rasterized to $domain grid (UAF-Rounce v.3), yearly',
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
  Code_Name	=> '$model\_$rcp\_c2_HiMAT_glaciers_area_y',
  Time_Series	=> 'yearly',
  Start_Date	=> '$dateBeg',
  End_Date	=> '$dateEnd',
  Name		=> 'HiMAT $model $rcp RGI-6.0 glacier area fraction rasterized to $domain grid (UAF-Rounce v.3), yearly',
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
');

#######################################################################

pp_done();
