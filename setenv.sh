#/bin/bash
#
# Set ALICE and shared libraries paths

export ALIBUILD_WORK_DIR=/home/user/alice/sw

alienv enter AliPhysics/latest-aliroot5-user

export LD_LIBRARY_PATH=./libs/RooUnfold/:${LD_LIBRARY_PATH}
