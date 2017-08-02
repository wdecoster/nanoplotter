# nanoplotter
This module provides functions for plotting data extracted from Oxford Nanopore sequencing reads and alignments, but some of it's functions can also be used for other applications.


## FUNCTIONS
* Check if a specified color is a valid matplotlib color `checkvalidColor(color)`  
* Check if a specified output format is valid `checkvalidFormat(format)`  
* Create a bivariate plot with dots, hexbins and/or kernel density estimates. Also arguments for specifying axis names, color and xlim/ylim. `scatter(x, y, names, path, color, format, plots, stat=None, log=False, minvalx=0, minvaly=0)`  
* Create cumulative yield plot and evaluate read length and quality over time `timePlots(df, path, color, format)`  
* Create length distribution histogram and density curve `lengthPlots(array, name, path, n50, color, format, log=False)`  
* Create flowcell physical layout in numpy array `makeLayout()`  
* Present the activity (number of reads) per channel on the flowcell as a heatmap `spatialHeatmap(array, title, path, color, format)`  


## INSTALLATION
```bash
pip install nanoplotter
```
or  
[![install with conda](https://anaconda.org/bioconda/nanoplotter/badges/installer/conda.svg)](https://anaconda.org/bioconda/nanoplotter)
```
conda install -c bioconda nanoplotter
```
