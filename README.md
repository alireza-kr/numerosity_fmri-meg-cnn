<a id="top"></a>
<div align="center">

<img src="https://github.com/alireza-kr/NuRBiM/raw/main/Files/Numerosity.jpg" alt="NuRBiM banner — a numerosity dot-array display beside a silhouette with an activated brain overlay" width="640">

*A numerosity display (left) and the brain activity it evokes — the two signals this project sets out to decode and compare.*

# 🧠 NuRBiM

### Numerosity Representation in Brain and Machines

How do brains — biological and artificial — represent *"how many"*?

[![License](https://img.shields.io/github/license/alireza-kr/NuRBiM?color=lightgrey)](LICENSE)
![MATLAB](https://img.shields.io/badge/MATLAB-Core%20Language-e16737)
![Python](https://img.shields.io/badge/Python-MNE--Python-3776AB?logo=python&logoColor=white)
[![Stars](https://img.shields.io/github/stars/alireza-kr/NuRBiM?style=social)](https://github.com/alireza-kr/NuRBiM/stargazers)
[![Forks](https://img.shields.io/github/forks/alireza-kr/NuRBiM?style=social)](https://github.com/alireza-kr/NuRBiM/forks)

🏷️ `cnn` `fmri` `meg` `numerosity`

**[Overview](#overview)** · **[Structure](#repository-structure)** · **[Usage](#how-to-use-this-repository)** · **[Data](#data-and-pretrained-models)** · **[Publications](#publications)** · **[License](#license)**

</div>

---

## 🧠 Overview <a id="overview"></a>

**NuRBiM** (Numerosity Representation in Brain and Machines) is the research codebase behind a PhD dissertation on one of the brain's most fundamental quantity senses: **numerosity** — the ability to judge "how many" items are in a scene at a glance, without counting.

This repository holds the full pipeline used to study that ability through three complementary lenses, from raw stimuli all the way to a brain-vs-machine comparison:

- 🧪 **Behavioral experiments** — psychophysics tasks built on dot-array stimuli
- 🧲 **Human neuroimaging** — parallel fMRI and MEG pipelines, from preprocessing to decoding
- 🤖 **Artificial neural networks** — the same representational framework applied to CNNs trained on numerosity

The work was carried out at the **Center for Mind/Brain Sciences (CIMeC), University of Trento**, and was awarded the 2026 [Glushko Dissertation Prize](https://cognitivesciencesociety.org/glushko-dissertation-prize/) by the Cognitive Science Society. 🏆

<br>

| | |
|---|---|
| 🔬 **Field** | Cognitive & computational neuroscience |
| 🧭 **Modalities** | fMRI · MEG · CNN (deep learning) |
| 💻 **Code** | 32 MATLAB scripts (~8,500 lines) + 3 Python notebooks |
| 🏆 **Award** | [Glushko Dissertation Prize](https://cognitivesciencesociety.org/glushko-dissertation-prize/) 2026, Cognitive Science Society |
| 🧾 **License** | CC0 1.0 — public domain |
| 🏛️ **Institution** | CIMeC, University of Trento |

---

## ✨ Highlights <a id="highlights"></a>

- 🧠 **Dual-modality neuroimaging** — parallel fMRI and MEG pipelines built on SPM, AFNI, FreeSurfer, and CoSMoMVPA
- 🔗 **MEG–fMRI fusion** — combines MEG's temporal resolution with fMRI's spatial resolution via representational fusion
- 📐 **Representational Similarity Analysis (RSA)** — model and neural representational dissimilarity matrices (RDMs), plus noise-ceiling estimation
- 🤖 **Brain-vs-machine comparison** — CORnet-Z convolutional networks trained on three numerosity datasets, with pretrained weights provided
- 🎯 **End-to-end experiment pipeline** — from dot-array stimulus generation to fMRI/MEG-ready trial sequences
- 👁️ **Eye-tracking from fMRI** — gaze estimation via DeepMReye to help control for oculomotor confounds

---

## 🗂️ Repository Structure <a id="repository-structure"></a>

```
NuRBiM/
├── Experiment/                Stimulus generation & experimental design (MATLAB)
│   ├── MakeExp.m                  Top-level script — builds the full stimulus set for a subject
│   ├── MakeDotSample.m            Generates "sample" dot-array stimuli
│   ├── MakeDotMatch.m             Generates "match" (comparison) stimuli
│   ├── MakeSeqFMRI.m              Builds fMRI trial / ISI sequences
│   ├── MakeSeqMEG.m               Builds MEG trial / ISI sequences
│   └── ...                        Geometry & sampling utilities
│
├── Analysis/                  fMRI, MEG, RSA & fusion pipeline
│   ├── RSA/
│   │   ├── MeasureModelRDM.m          Candidate model RDMs (DISTATIS, SVM-based, ...)
│   │   ├── MeasureNeuralRDM.m         Neural / DNN RDMs (multiple distance metrics)
│   │   └── MeasureNoiseCeiling.m      Noise-ceiling estimation for RSA
│   ├── mri_spm_*.m                fMRI preprocessing & 1st-level GLM (SPM)
│   ├── mri_afni_*.m               Surface mapping, smoothing & group stats (AFNI)
│   ├── mri_cosmo_*.m              ROI / searchlight decoding (CoSMoMVPA)
│   ├── mri_deepmreye.{m,ipynb}    Eye-tracking from fMRI (DeepMReye)
│   ├── meg_mne_*.ipynb            MEG preprocessing & source analysis (MNE-Python)
│   ├── meg_cosmo_*.m              Sensor-level GLM, decoding & time generalization
│   ├── meg_afni_mapSTC2SURF.m     MEG source-to-surface mapping
│   ├── meg_glm_source.m           Source-level GLM
│   └── fusion_cosmo.m             MEG–fMRI representational fusion
│
├── Files/
│   └── Numerosity.jpg         Banner figure used in this README
│
├── LICENSE                    CC0 1.0 Universal
└── README.md                  This file
```

<details>
<summary><strong>📋 Full script-by-script reference</strong> (click to expand)</summary>

#### `Experiment/`

| Script | Description |
|---|---|
| `MakeExp.m` | Top-level entry point — generates the complete stimulus set for one subject and saves it to a `.mat` file |
| `MakeDotSample.m` | Creates the "sample" dot-array stimulus for a trial |
| `MakeDotMatch.m` | Creates the "match" (larger / smaller comparison) stimulus for a trial |
| `MakeSeqFMRI.m` | Builds the stimulus and ISI (inter-stimulus interval) trial sequence for fMRI sessions |
| `MakeSeqMEG.m` | Builds the stimulus and ISI trial sequence for MEG sessions |
| `convVdPx.m` | Converts visual angle (degrees) to screen pixels |
| `convInx1d3d.m` | Converts between 1D and 3D array indices |
| `imcircle.m` | Draws the circular dots used in the numerosity arrays |
| `randDiamPick.m` | Randomly samples dot diameters under the experiment's constraints |
| `randPermPick.m` | Random permutation sampling helper |

#### `Analysis/` — fMRI (SPM · AFNI · FreeSurfer · CoSMoMVPA)

| Script | Description |
|---|---|
| `mri_spm_preprocessing.m` | Subject-level fMRI preprocessing (realignment, coregistration, segmentation, normalization, smoothing) |
| `mri_spm_1stLevel.m` | First-level (subject-level) GLM estimation |
| `mri_afni_makeSmooth.m` | Surface-based smoothing utilities |
| `mri_afni_makeSurface.m` | Cortical surface reconstruction / registration (AFNI + FreeSurfer) |
| `mri_afni_mapNII2SURF.m` | Projects volumetric statistical maps onto the cortical surface |
| `mri_afni_group.m` | Group-level analysis utilities |
| `mri_cosmo_glm_roi.m` | ROI-based GLM / pattern estimation |
| `mri_cosmo_glm_searchlight.m` | Whole-brain searchlight GLM / decoding |
| `mri_cosmo_group.m` | Group-level statistics for multivariate results |
| `mri_deepmreye.m` / `.ipynb` | Gaze / eye-position estimation directly from fMRI data via DeepMReye |

#### `Analysis/` — MEG (MNE-Python · AFNI · CoSMoMVPA)

| Script | Description |
|---|---|
| `meg_mne_preprocess.ipynb` | MEG preprocessing — filtering, artifact rejection, epoching |
| `meg_mne_source_analysis.ipynb` | MEG source reconstruction / localization |
| `meg_afni_mapSTC2SURF.m` | Maps MEG source estimates onto the cortical surface |
| `meg_cosmo_glm_sensor.m` | Sensor-level GLM analysis |
| `meg_cosmo_glm_sensor_searchlight.m` | Sensor-space searchlight GLM |
| `meg_cosmo_decoding_sensor.m` | Sensor-level multivariate decoding |
| `meg_cosmo_time_generalization_sensor.m` | Temporal generalization decoding across time |
| `meg_cosmo_vp_sensor.m` | Variance-partitioning analysis at the sensor level |
| `meg_cosmo_group.m` | Group-level statistics for MEG results |
| `meg_glm_source.m` | Source-level GLM analysis |

#### `Analysis/RSA/` — Representational Similarity Analysis

| Script | Description |
|---|---|
| `MeasureModelRDM.m` | Builds candidate model representational dissimilarity matrices (RDMs) |
| `MeasureNeuralRDM.m` | Computes neural / DNN RDMs under multiple distance metrics & cross-validation schemes |
| `MeasureNoiseCeiling.m` | Estimates the noise ceiling for RSA model comparisons |

#### Cross-Modal Fusion

| Script | Description |
|---|---|
| `fusion_cosmo.m` | Fuses MEG and fMRI representational geometries (individual-subject MEG × group-average fMRI) |

</details>

---

## 🧭 How to Use This Repository <a id="how-to-use-this-repository"></a>

This is research code accompanying a PhD thesis, not a plug-and-play package — there's no single command-line entry point. The pipeline is organized by modality and analysis stage:

1. **Generate stimuli** — run `Experiment/MakeExp.m` to build the dot-array stimulus set and the fMRI/MEG trial sequences for a subject.
2. **Preprocess** — `mri_spm_preprocessing.m` for fMRI, `meg_mne_preprocess.ipynb` for MEG.
3. **Model / decode** — subject-level GLMs and pattern estimation (`mri_spm_1stLevel.m`, `meg_cosmo_glm_sensor.m`, `mri_cosmo_glm_searchlight.m`, ...).
4. **Aggregate** — group-level statistics (`mri_cosmo_group.m`, `meg_cosmo_group.m`, `mri_afni_group.m`).
5. **Compare representations** — build RDMs and test them against model / DNN predictions (`Analysis/RSA/`).
6. **Fuse modalities** — combine MEG and fMRI representational geometries (`fusion_cosmo.m`).

> Most functions take a project-specific path/parameter structure (`mypath`, `mri`, `meg`, `mvpa`, ...) as input — define these for your own data layout before running. See [Toolboxes and Requirements](#toolboxes-and-requirements) below for each dependency's own setup instructions.

---

## 🧰 Toolboxes and Requirements <a id="toolboxes-and-requirements"></a>

| Toolbox | Used for |
|---|---|
| **MATLAB** | Core language for the experiment and analysis scripts |
| [**SPM**](https://www.fil.ion.ucl.ac.uk/spm/) | fMRI preprocessing & first-level GLM |
| [**AFNI**](https://afni.nimh.nih.gov/) | Surface mapping, smoothing & group-level analysis |
| [**FreeSurfer**](https://surfer.nmr.mgh.harvard.edu/) | Cortical surface reconstruction |
| [**CoSMoMVPA**](https://www.cosmomvpa.org/) | Multivariate pattern analysis — decoding, RSA, searchlight |
| [**Surfing**](https://surfing.sourceforge.net/) | Surface-based searchlight analysis |
| [**THINGSvision**](https://thingsvision.github.io/) | Extracting CNN / DNN activations |
| [**MNE-Python**](https://mne.tools/stable/index.html) | MEG preprocessing & source analysis |

The Python notebooks additionally rely on `numpy`, `pandas`, `scipy`, `scikit-learn`, and `nibabel`.

---

## 📦 Data and Pretrained Models <a id="data-and-pretrained-models"></a>

All data and model weights are hosted on [OSF](https://osf.io/) (Open Science Framework).

**Pretrained CORnet-Z weights**

| Model | Trained on | Link |
|---|---|---|
| CORnet-Z | DeWind dataset | [Download](https://osf.io/download/qdres/) |
| CORnet-Z | ISA2 dataset | [Download](https://osf.io/download/ek7vw/) |
| CORnet-Z | Natural dataset | [Download](https://osf.io/download/x748g/) |

**Stimulus image sets**

| Dataset | Link |
|---|---|
| DeWind | [Download (.zip)](https://files.osf.io/v1/resources/6gdfu/providers/osfstorage/66e16af872b893a38e459e9e/?zip=) |
| ISA2 | [Download (.zip)](https://files.osf.io/v1/resources/6gdfu/providers/osfstorage/66e16c8e5f653b21762e0442/?zip=) |
| Natural | [Download (.zip)](https://files.osf.io/v1/resources/6gdfu/providers/osfstorage/66e16cb2c372a19c5b1eb4ce/?zip=) |

---

## 📖 Publications <a id="publications"></a>

This repository accompanies the following work:

**PhD Dissertation**

> Karami, A. (2024). *The representation of numerosity in the human brain and machines* (Doctoral dissertation). Center for Mind/Brain Sciences (CIMeC), University of Trento. [doi.org/10.15168/11572_402591](https://doi.org/10.15168/11572_402591)

🏆 Winner of the 2026 [Glushko Dissertation Prize](https://cognitivesciencesociety.org/glushko-dissertation-prize/), Cognitive Science Society.

**Related Papers**

> Karami, A., Castaldi, E., Eger, E., & Piazza, M. (2025). Distinct neural representational geometries of numerosity in early visual and association regions across visual streams. *Communications Biology*, 8(1). [doi.org/10.1038/s42003-025-08395-z](https://doi.org/10.1038/s42003-025-08395-z)

> Karami, A., Truong, N., & Piazza, M. (2025). Investigation of Numerosity Representation in Convolution Neural Networks. *CCN 2025 Proceedings*. [doi.org/10.32470/4j01408](https://doi.org/10.32470/4j01408)

<details>
<summary><strong>📚 BibTeX</strong></summary>

```bibtex
@phdthesis{karami2024numerosity,
  title  = {The representation of numerosity in the human brain and machines},
  author = {Karami, Alireza},
  year   = {2024},
  school = {Center for Mind/Brain Sciences (CIMeC), University of Trento},
  doi    = {10.15168/11572_402591}
}

@article{karami2025geometries,
  title   = {Distinct neural representational geometries of numerosity in early visual and association regions across visual streams},
  author  = {Karami, A. and Castaldi, E. and Eger, E. and Piazza, M.},
  journal = {Communications Biology},
  volume  = {8},
  number  = {1},
  year    = {2025},
  doi     = {10.1038/s42003-025-08395-z}
}

@inproceedings{karami2025cnn,
  title     = {Investigation of Numerosity Representation in Convolution Neural Networks},
  author    = {Karami, A. and Truong, N. and Piazza, M.},
  booktitle = {CCN 2025 Proceedings},
  year      = {2025},
  doi       = {10.32470/4j01408}
}
```

</details>

---

## ⚖️ License <a id="license"></a>

Released under **CC0 1.0 Universal** — a public-domain dedication. You're free to use, modify, and redistribute this code with no restrictions. See [LICENSE](LICENSE) for the full text.

---

## 👤 Author and Acknowledgements <a id="author-and-acknowledgements"></a>

**Alireza Karami**<br>
PhD, Center for Mind/Brain Sciences (CIMeC), University of Trento

[![GitHub](https://img.shields.io/badge/GitHub-alireza--kr-181717?logo=github&logoColor=white)](https://github.com/alireza-kr)

This work builds on the open-source neuroimaging and machine-learning community — SPM, AFNI, FreeSurfer, CoSMoMVPA, Surfing, THINGSvision, and MNE-Python. 🙏

Questions or issues? Feel free to [open an issue](https://github.com/alireza-kr/NuRBiM/issues).

---

<div align="center">

**[⬆️ Back to top](#top)**

</div>
