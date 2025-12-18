#!/bin/csh

set baseDir=`pwd`

setenv HOME /sphenix/u/$LOGNAME
source /etc/csh.login
foreach i (/etc/profile.d/*.csh)
source $i
end
#source $HOME/.login
#source /opt/sphenix/core/bin/sphenix_setup.csh
source /sphenix/u/shanehorner/setup_sphenix_environment_csh
#printenv

# input parameter is process number
echo "input parameter: " $1

echo "running..."
cd ${baseDir}
pwd

set outfile = ${baseDir}/sim_rds/kshort_$1.root

rm -f $outfile

root -b -q  'roo_DeltaPAnalysis.C('$1')'


if (! -e $outfile || ! -s $outfile) then
    echo "ERROR: outfile missing or empty"
    exit 1
    endif


set size = `stat -c%s $outfile`
if ($size < 10000) then
    echo "file size too small - failing job. "
    exit 1
    endif

set outfile1 = ${baseDir}/sim_tgr/kshort_$1.root

rm -f $outfile1

root -b -q  'fast_bulk_plot_analyze_kshort_mom.C('$1')'


if (! -e $outfile1 || ! -s $outfile1) then
    echo "ERROR: plot outfile missing or empty"
    exit 1
    endif


set size1 = `stat -c%s $outfile1`
if ($size1 < 4000) then
    echo "plot file size too small - failing job. "
    exit 1
    endif
    
    




