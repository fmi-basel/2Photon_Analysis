# 2Photon_Analysis

# Introduction #
This code is analgous to the 1P extraction code, published alongside [this](https://www.science.org/doi/full/10.1126/science.abg7277) paper and on [GithHub](https://github.com/fmi-basel/1Photon_Analysis). It is also a wrapper that is build around [CNMF](https://www.sciencedirect.com/science/article/pii/S0896627315010843?via%3Dihub) and [NormCorre](https://www.sciencedirect.com/science/article/pii/S0165027017302753) and additionally features a SVM classification, which enables the automatic sorting of CNMF extracted components beyond the CNN included in the CNMF.

By now the CNMF package is also implemented in a workflow - [CaImAn](https://github.com/flatironinstitute/CaImAn), which also works with NormCorre and provides a maintained and updated package that is the state of the art for Ca-Extraction together with [Suite2P](https://github.com/MouseLand/suite2p). Therefore this code represents an adaptation that might help people to implement their own adaptations, but shouldn't replace usage of the well maintained packages already out there.

The code is also written to work only within a very specific Folder Structure, which is outlined in the pdf published alongside the code. It expects individual tif files or tif stacks that it concatenates and it works with single and multi-plane recordings.

# Installation #
The easiest way to use this code is to clone the repository in Matlab. You can check out how to do it [here](https://www.mathworks.com/help/matlab/matlab_prog/retrieve-from-git-repository.html) and [here](https://www.youtube.com/watch?v=O7A27uMduo0). Alternatively you can also download the code and save and maintain it locally on your computer. 

It requires a computer with a minimum of 64 GB RAM and for ease of processing should have at least 4 - 6 cores. It has been tested with Matlab 2017b and Matlab 2018a. It requires the following Matlab packages:  
  
* 'Optimization Toolbox'  
* 'Signal Processing Toolbox'  
* 'Image Processing Toolbox'  
* 'Statistics and Machine Learning Toolbox'  
* 'Curve Fitting Toolbox'  
* 'Parallel Computing Toolbox'  

To run the code, open the script *TwoP_Imaging_Analysis.m.m* in your editor and follow the instructions below.

# Usage #
## Parameter Selection ##
Parameters are described extensively in inline comments in the code and all set parameters should work in the specified range and can be changed if results do not match expectation, but are generally suited to obtain good first results. Moreover, as pointed out above, the code expects a certain folder structure that you need to adhere to.

## Session Selection and Motion Correction ##
After setting the parameters, running the code will open a prompt which asks you to add the top level folder of every Folder and the code will process every Session underneath the top level Folder and process them in order (Session_1 ... Session_n).
First, the data is converted to several .mat files, which are used for later processing stages. After this step, the code will pass the data to [NormCorre](https://github.com/flatironinstitute/NoRMCorre) and runs through all the files and generates a visualization for inspection - whih includes both a downsampled video and an enhanced projection image.
After finishing the Motion Correction all videos will have a downsampled video visualisation and max intensity projections for visual inspection. Here it is important that neurons in the center appear round with clear edges and that, focusing on major landmark in the video (i.e. blood vessels), there is no visible translational movement left. You can click through all video's by clicking enter after clicking into the command window. If you have accepted all sessions the code will automatically proceed to the next step. If you didnt't accept the motion correction, the code will ask you to draw a ROI in the enhanced correlation image, since most often the failure of motion correction can be explained by very bright spots in the FOV, which leads to a hang up of the MC. This step will repeat until the video is accepted or the maximum runtime is reached.

## CNMF ##
The code will then run all raw, motion corrected video's through an unedited version of [CNMF]([https://github.com/zhoupc/CNMF_E](https://github.com/flatironinstitute/CaImAn-MATLAB/tree/master)). Here an important feature is the selection of parallel pools that can speed up the process significantly, but also require large amounts of RAM, therefore a balance needs to be found. 

## CNMF - Postselection ##
During post-Selection the user will have the chance to sort all components that are not automatically excluded manually with a GUI that is slightly edited from the CNMFE implementation [CNMFE](https://github.com/zhoupc/CNMF_E). Follow the instruction displayed in the command window to proceed. If you set the SVM classifier variable to true, it will do the manual sorting step automatically and generate the final Output of each session.

All files will be saved alongside the raw data in a folder called processed_data_folder. To understand variables you can refer to the guide of CNMF, since variable naming is contained.

# References #
If you decide to use this code, please cite the relevant works especially of the people who wrote the original Code that runs this pipeline:    

Pnevmatikakis EA, Soudry D, Gao Y, Machado TA, Merel J, Pfau D, Reardon T, Mu Y, Lacefield C, Yang W, Ahrens M, Bruno R, Jessell TM, Peterka DS, Yuste R, Paninski L. [Simultaneous Denoising, Deconvolution, and Demixing of Calcium Imaging Data](https://www.sciencedirect.com/science/article/pii/S0896627315010843?via%3Dihub). Neuron. 2016 Jan 20;89(2):285-99. doi: 10.1016/j.neuron.2015.11.037. Epub 2016 Jan 7. PMID: 26774160; PMCID: PMC4881387.

Pnevmatikakis EA, Giovannucci A. [NoRMCorre: An online algorithm for piecewise rigid motion correction of calcium imaging data](https://www.sciencedirect.com/science/article/pii/S0165027017302753). J Neurosci Methods. 2017 Nov 1;291:83-94. doi: 10.1016/j.jneumeth.2017.07.031. Epub 2017 Aug 3. PMID: 28782629.

# Help #
If you have further questions, please reach out! You can find my contact information [here](https://julian-hinz.github.io/).


