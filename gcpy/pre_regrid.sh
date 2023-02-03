#!/bin/bash

# to combine the gridspec-create and esmf_regridweightgen into one step

srcdir=$(pwd -P)
mkdir -p ${srcdir}/grid_files
#mkdir -p ${srcdir}/grid_files/source_grid
#mkdir -p ${srcdir}/grid_files/destination_grid
# Define separator lines
thickline="\n===========================================================\n"
thinline="\n-----------------------------------------------------------\n"

cd ${srcdir}/grid_files
#-----------------------------------------------------------------
# Ask user to select source restart file grid type
#-----------------------------------------------------------------
printf "${thinline}Choose source restart file grid type:${thinline}"
printf "   1. gcs\n"
printf "   2. latlon\n"
printf "   3. sgcs\n"

valid_type=0
while [ "${valid_type}" -eq 0 ]; do
    read s_grid_type
    valid_type=1
    if [[ ${s_grid_type} = "1" ]]; then
	s_gridtype=gcs
    sweight=c
    elif [[ ${s_grid_type} = "2" ]]; then
    s_gridtype=latlon
    sweight=lat_lon_
    elif [[ ${s_grid_type} = "3" ]]; then
    s_gridtype=sgcs
    sweight=c
    else
        valid_type=0
	printf "Invalid grid type. Try again.\n"
    fi
done

#-----------------------------------------------------------------
# Ask user to select source restart file grid resolution
#-----------------------------------------------------------------
if [[ ${s_gridtype} = "gcs" ]] || [[ ${s_gridtype} = "sgcs" ]];then
    printf "${thinline}Enter source restart file grid resolution:${thinline}"
    read s_gridres
        if [[ ${s_gridtype} = "gcs" ]];then
        gridspec-create gcs ${s_gridres}
        else
        printf "${thinline}Enter the stretch factor${thinline}"
        read s_stretch_factor
        printf "${thinline}Enter the target latitude${thinline}"
        read s_target_lat
        printf "${thinline}Enter the target longitude${thinline}"
        read s_target_lon
        gridspec-create sgcs -s ${s_stretch_factor} -t ${s_target_lat} ${s_target_lon} ${s_gridres}
        fi

elif [[ ${s_gridtype} = "latlon" ]];then
    printf "${thinline}Choose source restart file grid resolution:${thinline}"
    printf "   1. 2x2.5\n"
    printf "   2. 4x5\n"

    valid_res=0
    while [ "${valid_res}" -eq 0 ]; do
        read s_grid_res
        valid_res=1
        if [[ ${s_grid_res} = "1" ]]; then
        s_gridres=2x2.5
    	gridspec-create latlon -b -180 -89.5 177.5 89.5 91 144
        sweightres=91x144
        elif [[ ${s_grid_res} = "2" ]]; then
        s_gridres=4x5
        gridspec-create latlon -b -180 -89 175 89 46 72
        sweightres=46x72
        else
            valid_res=0
    	printf "Invalid grid resolution. Try again.\n"
        fi
    done
fi

#cd ${srcdir}/grid_files/destination_grid
#-----------------------------------------------------------------
# Ask user to select destination restart file grid type
#-----------------------------------------------------------------
printf "${thinline}Choose destination restart file grid type:${thinline}"
printf "   1. gcs\n"
printf "   2. latlon\n"
printf "   3. sgcs\n"

valid_type=0
while [ "${valid_type}" -eq 0 ]; do
    read d_grid_type
    valid_type=1
    if [[ ${d_grid_type} = "1" ]]; then
	d_gridtype=gcs
    dweight=c
    elif [[ ${d_grid_type} = "2" ]]; then
    d_gridtype=latlon
    dweight=lat_lon_
    elif [[ ${d_grid_type} = "3" ]]; then
    d_gridtype=sgcs
    dweight=c
    else
        valid_type=0
	printf "Invalid grid type. Try again.\n"
    fi
done

#-----------------------------------------------------------------
# Ask user to select destination restart file grid resolution
#-----------------------------------------------------------------
if [[ ${d_gridtype} = "gcs" ]] || [[ ${d_gridtype} = "sgcs" ]];then
    printf "${thinline}Enter destination restart file grid resolution:${thinline}"
    read d_gridres

        if [[ ${d_gridtype} = "gcs" ]];then
        gridspec-create gcs ${d_gridres}
        else
        printf "${thinline}Enter the stretch factor${thinline}"
        read d_stretch_factor
        printf "${thinline}Enter the target latitude${thinline}"
        read d_target_lat
        printf "${thinline}Enter the target longitude${thinline}"
        read d_target_lon
        gridspec-create sgcs -s ${d_stretch_factor} -t ${d_target_lat} ${d_target_lon} ${d_gridres}
        fi

elif [[ ${d_gridtype} = "latlon" ]];then
    printf "${thinline}Choose source restart file grid resolution:${thinline}"
    printf "   1. 2x2.5\n"
    printf "   2. 4x5\n"

    valid_res=0
    while [ "${valid_res}" -eq 0 ]; do
        read d_grid_res
        valid_res=1
        if [[ ${d_grid_res} = "1" ]]; then
        d_gridres=2x2.5
    	gridspec-create latlon -b -180 -89.5 177.5 89.5 91 144
        dweightres=91x144
        elif [[ ${d_grid_res} = "2" ]]; then
        d_gridres=4x5
        gridspec-create latlon -b -180 -89 175 89 46 72
        dweightres=46x72
        else
            valid_res=0
    	printf "Invalid grid resolution. Try again.\n"
        fi
    done
fi

#cd ${srcdir}/grid_files

if [[ ${s_gridtype} = "gcs" || ${s_gridtype} = "sgcs" ]] && [[ ${d_gridtype} = "latlon" ]];then
    ESMF_RegridWeightGen -s ${sweight}${s_gridres}*gridspec* -d *${dweight}${dweightres}* -w ${sweight}${s_gridres}_to_${dweight}${d_gridres}.nc -m conserve
elif [[ ${s_gridtype} = "latlon" ]] && [[ ${d_gridtype} = "gcs" || ${d_gridtype} = "sgcs" ]];then
    ESMF_RegridWeightGen -s *${sweight}${sweightres}* -d ${dweight}${d_gridres}*gridspec* -w ${sweight}${s_gridres}_to_${dweight}${d_gridres}.nc -m conserve
elif [[ ${s_gridtype} = "gcs" || ${s_gridtype} = "sgcs" ]] && [[ ${d_gridtype} = "gcs" || ${d_gridtype} = "sgcs" ]];then
    ESMF_RegridWeightGen -s ${sweight}${s_gridres}*gridspec* -d ${dweight}${d_gridres}*gridspec* -w ${sweight}${s_gridres}_to_${dweight}${d_gridres}.nc -m conserve
else
    ESMF_RegridWeightGen -s *${sweight}${sweightres}* -d *${dweight}${dweightres}* -w ${sweight}${s_gridres}_to_${dweight}${d_gridres}.nc -m conserve
fi

#-----------------------------------------------------------------
# Done!
#-----------------------------------------------------------------
printf "\n${thinline}Generated!\n"
printf "\n  Weight file have been generated to /grid_files/${sweight}${s_gridres}_to_${dweight}${d_gridres}.nc${thinline}\n"

exit 0
    

