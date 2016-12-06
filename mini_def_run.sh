#!/bin/bash 
# This script attempts to run each of the versions in its own directory
# Using a decent choice of options (at this time)

# Change test behavior here:
# Revert back to programmed default variables by uncommenting the following line:
# USE_CLI_DEF=1
# Uncheck this for testing the ITAC binaries as well (if present)
TEST_ALSO_ITAC=1

EXECNAME="miniparsec" #how this benchmark is called right now
SUFFX="mpi"  #the suffix given to each exec in compile time
ITAC="itac"  #label for binaries compiled with ITAC support
FLAVOR=${1:-"bc bc.omp mpi2 mpi2.omp nbc nbc.omp"} #currently built flavors of the benchmark

PREFIX_COMMAND="mpirun" #you probably need this if you run on a PBS cluster

if [[ $USE_CLI_DEF -eq 1 ]]; then #integer comparison here folks
    CLILINE=""
else
    RMAX="65"       #the radius of the real-space sphere. Bigger=more ops
    GRIDSPACE="0.4" #spacing between grid points. Smaller=more ops
    VERBOSITY="1"   #how verbose should be the output? 0=no file, only std 1=normal 2=file for each mpi rank (PE)
    #NORDER="12"     #full order (front and back) of the stencil. Bigger=more data transfer, more ops.
    NORDER="10"     #full order (front and back) of the stencil. Bigger=more data transfer, more ops.
    CLILINE="$RMAX $GRIDSPACE $VERBOSITY $NORDER"
fi

#small archiving function if you are rerunning the test
function mini_arcfile () {
local copyfile=${1:-"miniparsec.out"}
if [[ -f $copyfile ]]
then
    local newname=$copyfile.$(date +%d%m%y-%H%M --reference $copyfile)
    if [[ -f $newname ]]
    then
        echo "mini_arcfile:: WARNING: $newname already exists, did not move"
    else
        cp -n $copyfile $newname
        echo "mini_arcfile:: copied $copyfile to $newname"
    fi
fi
}

RDIR=$(pwd)

if [[ $TEST_ALSO_ITAC -eq 1 ]]; then
    export VT_TIMER=gettimeofday #seems good enough 
#    export I_MPI_DEBUG=5
fi

for THISFLAVOR in $FLAVOR
do
    THISTEST="$EXECNAME.$THISFLAVOR.$SUFFX"
        if [[ -x "$THISTEST" ]]; then
            echo "found $THISTEST, will test now"
            mkdir -p t_$THISFLAVOR
            cd t_$THISFLAVOR
            mini_arcfile
            $PREFIX_COMMAND ../$THISTEST $CLILINE
            cd $RDIR
        fi

    if [[ $TEST_ALSO_ITAC -eq 1 ]]; then
        THISTEST="$EXECNAME.$ITAC.$THISFLAVOR.$SUFFX"
        if [[ -x "$THISTEST" ]]; then
            echo "found $THISTEST, will test now"
            mkdir -p t_$THISFLAVOR.$ITAC/trace-tmp/
            export VT_FLUSH_PREFIX=trace-tmp #we do not want to use /tmp for these
            cd t_$THISFLAVOR.$ITAC
            mini_arcfile
            $PREFIX_COMMAND ../$THISTEST $CLILINE
            cd $RDIR
        fi
    fi
done

echo "done with this script, go check the results on your own"
