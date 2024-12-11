#!/bin/bash
##
## psi_input_file_v1.7.1devel27b.sh
## 
## Made by Petar Marinkovic
## Login   <pmar@cake.lr.tudelft.nl>
## 
## Started on  Sun Jan 15 12:00:00 2005 Petar Marinkovic
## Last update Wed Jul 12 13:52:47 2006 Petar Marinkovic
##
## DESCRIPTION: Script that constructs input file fore
## ps_toolbox. Note that a threshold on delta_fDC is introduced
## input_project file, see input.template for more info.
##
## INPUT: [$1] : processing directory: location of interferograms
##        [$2] : output file
##        [$3] : doppler_threshold
##        [$4] : 'srp' or 'srd' (for h2ph)
##

## NOTE 1: this is development version, should be only used for
## testing!  The template is not compatible with previous releases of
## PSI scripts.
##
## NOTE 2: this is script is "crippled" version of psi_input_file.sh
## of fprits suite. All functions that would be sourced are hard coded
## within script
##
##
## CHANGE.LOG:
##
## [15.01.06, 12:00]: initial version
##
## [20.01.06, 12:00]: notation of storage conventions implemented, for
## more info see [data_storage_conventions.txt]
##
## [15.02.06, 10:37]: some profiling and input from Gini: (1) ifgs
## without h2ph do not break the loop, only warning now; (2)
## additional line for master_master.ifg
##
## [03.06.06, 11:22]: variables are grepped from res files with nawk,
##                    here defined as a function
##
## [03.07.06, 11:22]: orbit number listed in output.file
##
## [11.07.06, 11:08]: master_master entry dumped as last one
##
## [12.07.06, 12:43]: local version, no sourcing of functions nor
## dependence on fprits processing model
##
## [12.07.06, 13:33]: multiple master configuration possible
##
## [08.10.08, 13:33]: depsi v1.8.0 version
##
## TODO1: make code modular : organize some calls in a functions
##
## TODO2: general TODO for all the scripts, better processing overview

# $1 : input file : check
if [[ $# > 5 ]] || [[ $# < 4 ]]; then
    echo "Usage: $0 processing_directory supermaster output_file doppler_thresh srp/srd"
    echo ""
    echo "NOTE1: doppler_thresh is an optional argument; if no value" 
    echo "       specified all ifgs are accepted."
    echo ""
    echo "NOTE2: if there is more then one master in the stack " 
    echo "       output.file name will be extended with master orbit number."
    exit 127
fi

process=${1}
supermaster=${2}
srpd=${5}

if [ -d ${process}/sbas ]; then #for sbas processing
    process2=${process}/sbas/${supermaster}/rsmp/
    proc_flag=sbas #stb, redun processing
else
    process2=${process}
    proc_flag=sm #sm processing
fi

location=`pwd`
outputfile_orig=${location}/${3}
outputfile=${outputfile_orig}

# NOTE: 'new' outputfile is always created, also because of [>>]
if [ -f ${outputfile} ]; then
    rm -f ${outputfile}
fi

# doppler_threshold_parameter
# ---
doppler_thresh=${4}

# NOTE: doppler_thresh is an optional argument,it is read from
# input_parameters file, or if specified from command line as 2nd
# argument. Command line suppresses declaration in input_project. 
# If no value specified all ifgs are accepted.

if [ -z ${doppler_thresh} ]; then
    echo ""
    echo "WARNING: no doppler_thresh value specified"
    echo "WARNING: all ifgs will be accepted for psi"
    echo ""
fi

# - FUNCTIONS ------------- start -

# LIST ONLY DIRECTORY (requests (n)awk)
lsd ()
{
    # local argument
    ls -l ${1} | nawk '/^d/ {print $NF}'
}

# - FUNCTIONS -------------- stop -

# TODO: double loop in case of multiple masters specified

# get maters


#ifg_list=$(ls -d ${process}/*_*)

ifg_list=$(lsd ${process} | grep _ )

# counters initialization
fdc_rejcount=0   # 'rejection' on delfta_fDC counter
h2ph_rejcount=0  # counter for presence of h2ph

for ifg in ${ifg_list}
  do

  master=`echo ${ifg} | cut -f1 -d "_"`
  slave=`echo ${ifg} | cut -f2 -d "_"`

  ifg_dir=${ifg}
  res_file=${process}/${ifg_dir}/${ifg}.res
  
  if [ "${proc_flag}" = "sm" ]; then
      slave_res_file=${process}/${ifg_dir}/${slave}.res
  else
      slave_res_file=${process}/sbas/${supermaster}/rsmp/${ifg_dir}/${master}_rsmp.res 
      #yes, indeed master
  fi
  

  if [ ${master} -eq ${slave} ]; then
      
      date_master=`nawk '/First_pixel_azimuth_time/ {printf($3)}' \
          ${process}/${master}_${slave}/${master}.res`
      Bperp=0
      Btemp=0
      slcrsmp=${process2}${ifg_dir}/${slave}.rsmp 
      interf=dummy
      h2ph=dummy
      atmo=dummy
      doppler_master=0
      
      l0=`nawk '/First_line\ \(w.r.t.\ original_master\)/ {lines[last] = $0 } END { print lines[last] }' ${slave_res_file} | cut -f3`
      lN=`nawk '/Last_line\ \(w.r.t.\ original_master\)/ {lines[last] = $0 } END { print lines[last] }' ${slave_res_file} | cut -f3`
      p0=`nawk '/First_pixel\ \(w.r.t.\ original_master\)/ {lines[last] = $0 } END { print lines[last] }' ${slave_res_file} | cut -f3`
      pN=`nawk '/Last_pixel\ \(w.r.t.\ original_master\)/ {lines[last] = $0 } END { print lines[last] }' ${slave_res_file} | cut -f3`

      echo ""
      echo -e "IFG: ${master}_${slave} \t : OK (dummy)"
      
      # create temp variable to be parsed at the end
      #master_dummy=$(echo "${master} dummy ${date_master} dummy ${slcrsmp} ${interf} ${h2ph} ${atmo} ${Btemp} ${doppler_master}")
      master_dummy=$(echo "${master} ${date_master} ${slcrsmp} ${interf} ${h2ph} ${atmo} ${Bperp} ${Btemp} ${doppler_master} ${l0} ${lN} ${p0} ${pN}")
      # doppler_master, date_master added by FvL
      
      # check on H2PH output: DO NOT FORGET CROP ID
  elif [ -f ${process}/${ifg_dir}/*${srpd}.h2ph ] && [ ${master} -ne ${slave} ]; then
      
      date_master=`nawk '/First_pixel_azimuth_time/ {printf($3)}' \
          ${process}/${master}_${slave}/${master}.res`
      date_slave=`nawk '/First_pixel_azimuth_time/ {printf($3)}' \
          ${process}/${master}_${slave}/${slave}.res`
      
	# check on Dopplers, [NEW]: now using nawk
      
	# greps and converts to integers
      doppler_master=`nawk '/DC_const/ {printf("%d",$5)}' \
          ${process}/${master}_${slave}/${master}.res`
      
	# greps and converts to integers
      doppler_slave=`nawk '/DC_const/ {printf("%d",$5)}' \
          ${process}/${master}_${slave}/${slave}.res`
      
#	doppler_slave=$(grep DC_const \
#	    ${data}/${slave}/${crop}/${slave}.res \
#	    | sed 's/.*\ //; s/\ .*//')
#	
#	doppler_slave=${doppler_slave%%.*} # convert to integers, easier
#	
#	echo ${doppler_slave}       # visual.debug
#	echo ${doppler_master}      # visual.debug
	
	# NOTE1: expr doesn't do floating point, that's why bc
	# NOTE2: is conversion to integers now necessary?? NO!
	# NOTE3: scale apply only to division and extra space in case of --
      
      doppler_diff=$(echo "scale=2 ; ${doppler_master} - ${doppler_slave}" | bc)
      doppler_diff_orig=${doppler_diff} #added by FvL
      
	# echo ${doppler_diff}        # visual.debug 
      
	# make doppler_diff absolute value or in 2nd [if] of loop
	# double check, abs value is much easier
      if [[ ${doppler_diff} -lt "0" ]]; then
	  doppler_diff=$(echo "scale=2 ; ${doppler_diff} * -1" | bc)
      fi

        # Note: &&, ||, <, and >  operators works within [[ ]] test
      if [[ ${doppler_diff} -lt ${doppler_thresh} ]]; then
	  
	    # again grep with awks: more robust since doris > v3.16
	    # has a bit different res file structure
	  
	  Bperp=`nawk '/perp/ {printf( $3 )}' ${res_file}`

	  Btemp=`nawk '/temp/ {printf( $3 )}' ${res_file}`
	    # Btemp=`grep temp ${res_file} | sed 's/.*\ //; s/\ .*//'`
	    # and convert to years
	  Btemp=$(echo "scale=9 ; ${Btemp} / 365" | bc)
	  
	  if [ "${proc_flag}" = "sm" ]; then
	      slcrsmp=${process2}${ifg_dir}/${slave}.rsmp 
	      
	      interf=${process2}${ifg_dir}/${master}_${slave}.${srpd}
	      
	      h2ph=${process2}${ifg_dir}/${master}_${slave}_${srpd}.h2ph
	      
	      atmo=${process2}${ifg_dir}/${master}_${slave} # dummy
	  else
	      slcrsmp=${process2}${ifg_dir}/${ifg}.${supermaster}.rsmp 
	      
	      interf=${process2}${ifg_dir}/${ifg}.${supermaster}.${srpd}
	      
	      h2ph=${process2}${ifg_dir}/${ifg}.${supermaster}_${srpd}.r4.h2ph
	      
	      atmo=${process2}${ifg_dir}/${ifg}.${supermaster} # dummy
	  fi

          l0=`nawk '/First_line\ \(w.r.t.\ original_master\)/ {lines[last] = $0 } END { print lines[last] }' ${slave_res_file} | cut -f3`
          lN=`nawk '/Last_line\ \(w.r.t.\ original_master\)/ {lines[last] = $0 } END { print lines[last] }' ${slave_res_file} | cut -f3`
          p0=`nawk '/First_pixel\ \(w.r.t.\ original_master\)/ {lines[last] = $0 } END { print lines[last] }' ${slave_res_file} | cut -f3`
          pN=`nawk '/Last_pixel\ \(w.r.t.\ original_master\)/ {lines[last] = $0 } END { print lines[last] }' ${slave_res_file} | cut -f3`

	  echo ""
	  echo -e "IFG: ${master}-${slave} \t : OK : \t [\t Doppler.diff: ${doppler_diff} \t]"
	  
            # ----------------------
            # DUMP to log/outputfile
	  #echo ${master} ${slave} ${date_master} ${date_slave} ${slcrsmp} ${interf} ${h2ph} ${atmo} ${Btemp} ${doppler_diff_orig} ${l0} ${lN} ${p0} ${pN}
	  echo ${slave} ${date_slave} ${slcrsmp} ${interf} ${h2ph} ${atmo} ${Bperp} ${Btemp} ${doppler_diff_orig} ${l0} ${lN} ${p0} ${pN} >> ${outputfile}

      else
	  
	  echo ""
	  echo -e "IFG: ${master}-${slave} \t : REJECTED (dop)"
	  
	  fdc_rejcount=$(( ${fdc_rejcount} + 1 ))
	  
	  continue
	  
      fi # doppler check and output
      
  else
      
      echo ""
      echo -e "IFG: ${master}-${slave} \t : REJECTED (h2ph)"
      
      h2ph_rejcount=$(( ${h2ph_rejcount} + 1 ))
      
      continue
      
  fi # check on doppler and h2ph
  
  
# NOTE: since master_master entry needs to be last one, dumping moved
# to main loop/check
#    # ---------------------- 
#    # DUMP to log/outputfile
#
#    echo ${orbit} ${slcrsmp} ${interf} ${h2ph} ${atmo} ${Btemp} \
#	>> ${log}/${outputfile}

done # loop on ifg

# ------------------------------------------
# NOTE: master_master entry MUST BE LAST ONE!
  
  echo ${master_dummy} >> ${outputfile} # variable created in main loop

# ------------------------------------------

  # Loop summary for ${master}
  
  echo ""
  echo "-------------------------------------"
  echo -e "Summary for master: \t [${supermaster}]"
  echo "-------------------------------------"
  echo ""

  echo "---"
  echo "INFO:"
  echo ""
  echo -e "Total number of accepted ifgs: \t $(( `lsd ${process} | grep _ | cut -f2 -d "_" | wc -l` - ${fdc_rejcount} - ${h2ph_rejcount})) \t [including master-master]"
  echo -e "Total number of rejected ifgs: \t $(( ${fdc_rejcount} + ${h2ph_rejcount} ))"
  echo ""
  echo -e "Name of input/output file: \t `basename ${outputfile}`"
  echo -e "In dir: \t \t \t `pwd`"
  echo "---"
  echo ""

  if [ ${fdc_rejcount} -ne 0 ]; then
      echo "---"
      echo "WARNING: Some ifgs failed because of fDC thresh"
      echo "-"
      echo -e "Delfta_fDC thresh: \t ${doppler_thresh}"
      echo -e "No of ifgs rejected: \t ${fdc_rejcount}"
      echo "---"
      echo ""
  fi

  if [ ${h2ph_rejcount} -ne 0 ]; then
      echo "---"
      echo "WARNING: Some ifgs failed because of no h2ph file"
      echo "-"
      echo -e "No of ifgs rejected: \t ${h2ph_rejcount}"
      echo "---"
      echo ""
#      echo "How did you compute reference phase? Input card: EXACT?"
#      echo "Was ifg computed, check coregistration?"
#      echo ""
  fi
  
exit 0

# EOF
