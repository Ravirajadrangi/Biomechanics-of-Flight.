# ME303 - Win 17 Final Project
_(Updated: 03/13/17)_

## Perspectives
This repository stores the analysis codes and physical backgrounds of my final project for the Stanford ME303 course (Winter 2017): Biomechanics of Flight. All material can be freely distributed with proper citation of this repository: [https://github.com/Mipanox/ME303](https://github.com/Mipanox/ME303). 

## Introduction
The project aims at exploring birds' behavior in landing, in particular the decelaration phase of an american avocet and a northern shoveler. In the [videos](https://github.com/Mipanox/ME303/tree/master/kinematics/videos) (or on Nimia: [avocet](https://app.nimia.com/video/712475/american-avocets-take-off-and-landing/); [shoveler](https://app.nimia.com/video/723940/northern-shoveler-water-skiing/)), the two birds exhibit quite distinct ways of braking: The avocet uses its wings as a "parachute" whereas the shoveler water-skies and in the same time flaps the wings back-and-forth. Simple models are applied to describe the physics behind the behavior. 

_(Note: As of this writing the author has not included any conclusion remarks in any of the notebooks. Will be there soon)_

## Kinematics
The extracted kinematics of representatitve parts of the body motions (e.g. wings, tails, torso) is shown [here](https://nbviewer.jupyter.org/github/Mipanox/ME303/blob/master/kinematics/kinematics.ipynb). This notebook includes a brief overview of the videos, some crude calibrations, and then presents the smoothed kinematic curves: for instance the bulk velocity variation with time.

## Dynamics
### Assumptions
In this [notebook](https://nbviewer.jupyter.org/github/Mipanox/ME303/blob/master/dynamics/models.ipynb) the simplified models adopted in the simulation ale explained. Some qualitative descriptions of the assumptions are made, but refer to the [estimates](https://nbviewer.jupyter.org/github/Mipanox/ME303/blob/master/others/estimates.ipynb#Assumption-Justification) for justification (or not).

### Parameters
From the kinematics, it is possible to retrieve useful quantities to be thrown into the models. This is illustrated [here](https://nbviewer.jupyter.org/github/Mipanox/ME303/blob/master/others/estimates.ipynb).

### Models - Simulations
The whole process of landing is modeled by three mechanisms: (1) [Steady flight](https://nbviewer.jupyter.org/github/Mipanox/ME303/blob/master/dynamics/QS.ipynb), (2) [parachuting](https://nbviewer.jupyter.org/github/Mipanox/ME303/blob/master/dynamics/Brake.ipynb), and (3) [back-and-forth flapping plus water-skiing](https://nbviewer.jupyter.org/github/Mipanox/ME303/blob/master/dynamics/Planing.ipynb).
Based on the assumptions, the birds are fully characterized by point masses, (half-)elliptical wings, and flat plates for the feet, with certain properties inferred from the kinematics/biology. The [codes](https://nbviewer.jupyter.org/github/Mipanox/ME303/tree/master/codes/) simulate the forces and powers that would be generated for the modelled "artificial" birds. Comparisons with observations are also conducted and examined in each notebook.

## Conclusion
_Coming soon..._