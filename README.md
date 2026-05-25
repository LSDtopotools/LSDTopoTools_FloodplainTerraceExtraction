# LSDTopoTools_TerraceExtraction

**THIS REPOSITORY IS RETIRED AND FOR ARCHIVE PURPOSES ONLY**, For loatest versions see [lsdtopotools2](https://github.com/LSDtopotools/LSDTopoTools2). *lsdtt3*, the latest version is under testing and is due to be released in Autumn 2026. 

This repository contains code for extracting terraces automatically from DEMs using slope and channel relief thresholds.  The methodology is outlined in [Clubb et al. (2017)](https://doi.org/10.5194/esurf-2017-21).

The user provides latitude and longitude coordinates of two points along a channel of interest.  The code then extracts a swath profile along the channel between these two points following [Hergarten et al. (2014)](http://www.earth-surf-dynam.net/2/97/2014/esurf-2-97-2014.html).  The terraces are extracted from the swath profile based on thresholds of local gradient and elevation compared to the baseline channel.  These thresholds are calculated statistically from the DEM using quantile-quantile plots, and do not need to be set by the user for the landscape in question.

For more details on running the code please see the [LSDTopoTools documentation](https://lsdtopotools.github.io/LSDTT_documentation/LSDTT_floodplains_terraces.html)
