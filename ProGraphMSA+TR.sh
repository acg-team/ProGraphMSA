#!/bin/bash
BINARY="bin/ProGraphMSA_32"
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

if [[ ! -x "${BINARY}" ]]; then
   echo "ProGraphMSA binary not found at \"${BINARY}\". Make sure you run the program from the installation directory"
   exit 1
fi

if [[ ! $* =~ trust2treks ]]; then
   # Check whether T-REKS is installed and download if necessary
   if [[ ! -f "T-Reks.jar" ]]; then
      read -p "T-REKS does not seem to be installed here. Download it? [y/n]" -n 1 INST
      echo
      if [[ ${INST} =~ ^[Yy]$ ]]; then
         wget http://bioinfo.montp.cnrs.fr/t-reks/T-Reks.jar
      fi
   fi
else
   # Check whether TRUST is installed and download if necessary
   if [[ ! -f "Align/nl/vu/cs/align/SelfSimilarity.class" ]]; then
      read -p "TRUST does not seem to be installed here. Download it? [y/n]" -n 1 INST
      echo
      if [[ ${INST} =~ ^[Yy]$ ]]; then
         wget http://www.ibi.vu.nl/programs/trustwww/trust.tgz
         tar xzf trust.tgz
      fi
   fi
fi

${BINARY} --repeat_indel_rate 0.1 --repeat_indel_ext 0.3 --mldist --repeats --fasta $@
#--custom_tr_cmd trust2treks.py
#
