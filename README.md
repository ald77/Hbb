Hbb
===

Shared code for Hbb SUSY analysis.

Code may be compiled with ./compile.sh. Executables and scripts are stored in scripts directory and are all intended to be run from the root Hbb directory (i.e. not from within the scripts directory).

Repository is setup with the assumption that (assuming the repository is locally called Hbb) there is a directory Hbb/../data. This directory is intended to store cfA-style .root files containing data from cfA needed in the analysis. The directory can be filled (but not created) by running ./scripts/make_skims.sh, ./scripts/fix_skims.sh, and then ./scripts/shrink_skims.sh. This process will be simplified soon.
