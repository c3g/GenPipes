PCGR_VERSION="1.2.0"
OUTPUT_DIRECTORY="$MUGQIC_INSTALL_HOME/software/pcgr/tmp/PCGR"

git clone \
  -b "v${PCGR_VERSION}" \
  --depth 1 \
  https://github.com/sigven/pcgr.git \
  "${OUTPUT_DIRECTORY}"
