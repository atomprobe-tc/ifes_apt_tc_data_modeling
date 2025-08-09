# NeXus NXapm

## Towards a global data model for atom probe research

Recently, coordinated efforts that have build on the collective knowledge of the atom probe community have resulted
in the development of a covering data model for atom probe tomography and related field ion microscopy. This model
uses the [NeXus data modeling framework](https://www.nexusformat.org). The proposed data model `NXapm` captures
all aspects of data acquisition and data analysis to arrive at calibrated reconstructed datasets with applied
ranging definitions with transparent communication about peak fitting and analysis routines.
The proposal has recently been proposed for standardization with the [NeXus international Advisory Committee (NIAC)](https://www.nexusformat.org/NIAC.html).
The proposal was successful resulting in `NXapm` since July 2025 being an official part of the NeXus standard.
[This is the respective pull request.](https://github.com/nexusformat/definitions/pull/1422)

The activities have been acknowledged by key members of the atom probe community thus qualifying as suggested
to be used [global reporting standard for atom probe.](https://doi.org/10.1093/mam/ozae081)

With [pynxtools-apm](https://www.github.com/FAIRmat-NFDI/pynxtools-apm.git) a reference implementation exists that maps file formats in atom probe except for the CamecaROOT-based ones to `NXapm`.

