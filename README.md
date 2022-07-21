## Model-based hearing-enhancement strategies for cochlear synaptopathy pathologies

**Drakopoulos, F., Vasilkov, V., Osses Vecchi, A., Wartenberg, T. & Verhulst, S. Model-based hearing-enhancement strategies for cochlear synaptopathy pathologies. Hearing Research 108569 (2022). https://doi.org/10.1016/j.heares.2022.108569**

This repository contains the presented hearing-loss compensation strategies for cochlear synaptopathy (CS), as well as the evaluation results of our participants and the scripts that were used to record and analyse EEGs. The supporting paper can be found [here](https://doi.org/10.1016/j.heares.2022.108569) and a pre-print version is available under the *doc* folder (https://www.biorxiv.org/content/10.1101/2022.01.10.475652v2). The hearing-enhancement algorithms were designed to compensate for the degraded AN responses of CS-affected auditory peripheries, based on simulated outcomes of a biophysically inspired auditory model ([Verhulst et al. 2018, v1.2](https://github.com/HearingTechnology/Verhulstetal2018Model)). 

### Hearing-enhancement algorithms

The *algorithms* folder includes the implementations of the hearing-enhancement algorithms in MATLAB. Three processing functions are included (`g_70dB.m`, `gm_70dB.m`, and `gmref_70dB.m`), which correspond to the implementations of the original processing strategy (Eq. 2), the modified processing strategy (Eq. 9), and the clean-envelope (reference) modified strategy, respectively, as presented in the paper. Each strategy can be applied to an input signal to improve temporal-envelope processing for 3 CS profiles: *13,0,0* (loss of LSR and MSR ANFs), *10,0,0* (loss of LSR and MSR ANFs, 23% loss of HSR ANFs) and *7,0,0* (loss of LSR and MSR ANFs, 46% loss of HSR ANFs). 

### Experimental evaluation

The developed algorithms were evaluated in normal-hearing (NH) subjects of two age groups: Young NH (yNH) and older NH (oNH) subjects, without and with suspected age-related CS, respectively. The *data* folder contains the evaluation results of our participants, which include measured audiometric thresholds, envelope-following responses (EFRs), amplitude-modulation (AM) detection thresholds, and speech intelligibility in terms of speech-recognition thresholds (SRTs) and word-recognition scores (WRSs).

### EEG recording and analysis 

The *EEG_recording* folder contains all the necessary files to record EFRs on our BioSemi EEG setup. The `run_experiment.m` script is used to record EEG responses to a specific stimulus, implemented in MATLAB using [Playrec](https://github.com/PlayrecForMatlab/playrec). The two stimuli that were used for our experimental evaluation (SAM tone and HP-filtered speech) are provided under the *stimuli* folder. In the `run_experiment.m` script, the playback parameters (*pars*) need to be defined depending on the stimulus, given as an argument or adapted inside the script by uncommenting the relative sections. To generate a processed stimulus using one of the 3 CS-compensating algorithms, the corresponding number of the desired CS profile (*1300*, *1000* or *700*) can be given as the third argument. 

Additionally, the *EFR_analysis* folder contains the `analyse_bdf.m` script (along with all the supplementary functions) that was used to extract and analyse the EFRs from the recorded bdf files of the Biosemi setup.

----
## Citation
If you use this code, please cite the corresponding paper:

Drakopoulos, F., Vasilkov, V., Osses Vecchi, A., Wartenberg, T. & Verhulst, S. Model-based hearing-enhancement strategies for cochlear synaptopathy pathologies. Hearing Research 108569 (2022). https://doi.org/10.1016/j.heares.2022.108569

##
For questions, please reach out to one of the corresponding authors:

* Fotios Drakopoulos: fotios.drakopoulos@ugent.be
* Sarah Verhulst: s.verhulst@ugent.be

> This work was funded with support from the EU Horizon 2020 programme under grant agreement No 678120 (RobSpear).

