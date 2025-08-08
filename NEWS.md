# ggdmc News

## Version 0.2.8.9 (Release Date)

### New Features
- Added an external interface to automatically start three (or more) independent sampling processes.
- Reorganise ggdmc to ggdmcHeaders, ggdmcPrior, ggdmcLikelihood, and ggdmcModel.
- Separate specific choice response-time models to their package and link their likelihood to the ggdmc with the ggdmcHeaders.
- Replace the ggplot2 plotting code with lattice package to make the installation easier.

### Bug Fixes
- Fixed an issue where `migration` sampler misses one chain when only a small number of 
chains was sampled for migration. 
- Resolved a bug in `other_function()` causing [describe problem]. (#PR, @contributor)



