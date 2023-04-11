#setting home alias for x17 directory
alias home='cd /storage/projects/utefx17/martin;pwd'

#garfield directory
export GARFIELD_HOME='/storage/projects/utefx17/SourceCode/garfield'

#source adding modules
source /etc/profile.d/20_meta_modules.sh

#adding cmake
#module add cmake-3.14.5
module add cmake/cmake-3.17.3-gcc-8.3.0-z6akqlo

#sourcing Geant4
source /storage/projects/utefx17/SourceCode/geant4/geant4-v11.0.3-install/bin/geant4.sh

#sourcing Garfield
source $GARFIELD_HOME/install/share/Garfield/setupGarfield.sh

#for qt building
export LLVM_INSTALL_DIR=/storage/projects/utefx17/SourceCode/dependencies/llvm-14/

# If not running interactively, don't do anything
case $- in
    *i*) ;;
      *) return;;
esac

#sourcing ROOT (causes annoying version info errors)
source /storage/projects/utefx17/SourceCode/ROOT/install/bin/thisroot.sh

#trimming long directory output
export PROMPT_DIRTRIM=1