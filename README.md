# Calculate Size of Nucleosome Free Regions (NFR)
This command line tool provides an easy interface for calculating the size of NFRs.
It is based on a simple algorithm that can be summarised as follows:

1) Take the TSS as a reference position. They are given by the passed `.bed` or `.gtf` file 
2) Consider a range around the TSS (`x_range`), fetch MNase data values, and smooth it. 
3) Calculate derivative and determine minima and maxima positions.
4) The position of the NFR is defined as the lowest local minima in the upstream region
5) The +1 and -1 nucleosome positions are subsequently defined by the maxima that sandwich the NFR
6) If there is a peak that is significantly larger (defined by  a threshold), then the nucleosome position is accordingly replaced

The detection is considerably robust and detects the NTF even when the given TSS is not arount the actual +1 nucleosome.

| |  | |
:-------------------------:|:-------------------------:|:-------------------------:
![example 1](figures/arr1.png)|![example 2](figures/ypr202w.png) |![example 3](figures/ypr204w.png)  


## Requirements
The program requires Python3 and pip to be installed

## Installation
Install the dependencies through

```commandline
python3 -m pip install -r requirements.txt
```

## Usage
The MNase signal should be converted to a bigwig `.bw` file. Gene annotations can be passed as a `.bed` por `.gtf` file.
NFR sizes are saved as a `.tsv` file together with the corresponding annotation names. 

```commandline
python3 nfr.py --bw=path/to/bw/file --annot=path/to/annotation [--smooth_ws=50 --sig_nfold=3 --mind_nfr=50 --maxd_total=500 --out_path=path/to/output/file --verbosity=0]
```

Parameters in brackets are optional.
- `--smooth_ws`: window size of smoothing window that is used for filtering the signal to decrease the impact of noise.
- `--sig_nfold`: nfold difference to initially computed position to consider another peak as actual +1/-1 nucleosome
- `--mind_nfr`: minimum difference from TSS to NFR
- `--maxd_total`: maximum distance from the TSS that is considered for finding +1/-1 nucleosome. Thus, the total window size is 2 * `maxd_total`
- `--out_path`: output directory
- `--verbosity`: verbosity flag to set amount of visualisation during execution. 