# Dregon
Official repository of the SPCUP19: DREGON CHALLENGE

## Usage
1. Download or clone this repository
2.  Download the development data ('dev\_flight.mat' and 'dev\_static.mat') from the [resources page](https://piazza.com/ieee_sps/other/spcup2019/resources) on the PIAZZA platform.
3. Unzip the content of the zip files into the code's folder (see next section)
4. Run the main file 'baseline.m' with MATLAB.

__Note__: baseline.m features hard-coded path for UNIX-like system. For Windows system, a modification of this variables in needed.
__Note__: the baseline code has been tested on MATLAB R2017a.

## Files organization
After 3. point of the previous section, the folder should be organized as follows:
* SPCUP19/
    * MBSSLocate/
        * ...
        * ... files and folders of [MBSS Locate](http://bass-db.gforge.inria.fr/bss_locate/) MATLAB toolbox
    * `baseline.m`
    * dev_flight/
        * audio/
        * `SPCUP19\_dev\_flight.mat`
    * dev_static/
        * audio/
        * `SPCUP19\_dev\_static.mat`

## Baseline with MBSS Locate
As baseline MBSS Locate is used. It is an implementation of the state-of-the-art _steered response power with phase transform_ (`SRP-PHAT`) algorithm (Dibiase et al., 2001). This implementation is freely available, together with 7 other _angular spectrum-based_ localization techniques (Blandin et al., 2012), in a Matlab toolbox named Multichannel BSS Locate  [here](http://bass-db.gforge.inria.fr/bss_locate/).
In this challenge, the _generalized cross-correlation with phase transform_ (`GCC-PHAT`) is used for estimating the the angular spectrum of each pair of microphones.

The version shipped within this repository is the base one. A full version with example and data is available on the toolbox website as well.

This toolbox is easy to use:
1. a few line of code in the main files (from 107 till 140) set the parameters for the signals (e.g. sampling frequency, frame/block size), the source (e.g. static/moving, single/multiple) the microphones array (e.g. the array geometry and moving/static behaviour) and the localization method (e.g. angular spectrum technique)
2. the function `MBSS\_InputParam2Struct(...)` tokes as input the parameters obove mentioned and return a structure which is passed to the main fuction of the toolbox
3. `MBSS\_locate\_spec(...)` performs the localization, that is the estimation of both the azimuth(s) and elevation(s) of the target source(s) from the multichannel signals.

__Note__: MBSS Locate require the knownledge of the microphone array.

## References
* DiBiase, J. H., _A High Accuracy, Low-Latency Technique for Talker Localization in Reverberant Environments using Microphone Arrays_, (Ph.D.). Brown Univ, 2000,
* C. Blandin, A. Ozerov and E. Vincent, _Multi-source TDOA estimation in reverberant audio using angular spectra and clustering_, Signal Processing 92, pp. 1950-1960, 2012.
