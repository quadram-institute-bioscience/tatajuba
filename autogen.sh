#!/bin/sh

# developers can point to a common location, not as submodule but end users, if installing properly with git clone --recursive, 
# need to point to submodule
if [ ! -e "./biomcmc-lib" ]; then  # -e is file (-f), directory (-d), link (-h) etc.
  ln -s submodules/biomcmc-lib ./biomcmc-lib
fi

export AUTOMAKE="automake --foreign -a"
autoreconf -f -i
