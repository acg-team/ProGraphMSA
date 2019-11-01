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

# Check whether ProGraphMSA is run from the correct directory

if [[ ! -x "${BINARY}" ]]; then
   echo "ProGraphMSA binary not found at \"${BINARY}\". Make sure you run the program from the installation directory"
   exit 1
fi

${BINARY} --darwin --cs_profile 3rd_party/K4000.lib --mldist --estimate_aafreqs --fasta $@
