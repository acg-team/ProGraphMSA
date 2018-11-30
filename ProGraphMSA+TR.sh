#!/bin/bash

# Get the base directory of the argument.
# Can resolve single symlinks if readlink is installed
function scriptdir {
    cd "$(dirname "$1")"
    cd "$(dirname "$(readlink "$1" 2>/dev/null || basename "$1" )")"
    pwd
}
# Directory of this script
DIR="$(scriptdir "$0" )"

BINARY="ProGraphMSA"
SYSTEM=$(uname -s)
ARCH=$(uname -m)

if [[ "$SYSTEM" = "Linux" ]]; then
   if [[ "$ARCH" = "i386" ]]; then
      BINARY="bin/ProGraphMSA_32"
   elif [[ "$ARCH" = "x86_64" ]]; then
      BINARY="bin/ProGraphMSA_64"
   fi
elif [[ "$SYSTEM" = "Darwin" ]]; then
   BINARY="bin/ProGraphMSA_osx"
fi

if [[ "$PROGRAPHMSA_PREFIX" != "" ]]; then
   BINARY="${PROGRAPHMSA_PREFIX}/${BINARY}"
elif [[ -x "BUILD/src/ProGraphMSA" ]]; then
   BINARY="BUILD/src/ProGraphMSA"
fi

# Check whether ProGraphMSA+TR is run from the correct directory

if [[ ! -x "${DIR}/${BINARY}" ]]; then
   echo "ProGraphMSA binary not found at \"${DIR}/${BINARY}\". Make sure you run the program from the installation directory"
   exit 1
fi

if [[ ! $* =~ trust2treks ]]; then
   # Check whether T-REKS is installed and download if necessary
   if [[ ! -s "$DIR/T-Reks.jar" ]]; then
      read -p "T-REKS does not seem to be installed. Accept non-commercial license and download it? [y/n]" -n 1 INST
      echo
      if [[ ${INST} =~ ^[Yy]$ ]]; then
         # Old location: http://bioinfo.montp.cnrs.fr/t-reks/T-Reks.jar
         # Might have to use --no-check-certificate or install the proper CA for the current CNRS
         wget -O "$DIR/T-Reks.jar" https://dali.crbm.cnrs.fr/tools/treks/T-Reks.jar || exit
      fi
   fi
else
   # Check whether TRUST is installed and download if necessary
   if [[ ! -s "$DIR/Align/nl/vu/cs/align/SelfSimilarity.class" ]]; then
      read -p "TRUST does not seem to be installed here. Accept TRUST license and download it? [y/n]" -n 1 INST
      echo
      if [[ ${INST} =~ ^[Yy]$ ]]; then
         wget -O "$DIR/trust.tgz" http://www.ibi.vu.nl/programs/trustwww/trust.tgz || exit
         tar xzf "$DIR/trust.tgz" -C "$DIR" || exit
      fi
   fi
fi

exec "${DIR}/${BINARY}" --repeat_indel_rate 0.1 --repeat_indel_ext 0.3 --mldist --repeats --fasta $@
#--custom_tr_cmd trust2treks.py
#
