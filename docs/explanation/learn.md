# Exchange

## Calling for more exchange about (meta)data in the field

The proprietary software package IVAS now APSuite by AMETEK/Cameca is the workhorse
for data acquisition and analysis in the atom probe community. Several file formats
that this software uses are proprietary with the `.apt` file format as one exception.
The idea of the `.apt` file format is opening up atom probe data from Cameca instrument,
i.e., Local Electrode Atom Probe (LEAP) instruments for scientific data analysis and
public dissemination. The International Field Emission Society (IFES) and Cameca have
worked together to communicate a documentation for the format which enabled the
community to develop open-source reading capabilities as implemented also in the
`ifes_apt_tc_data_modeling` library through its `apt` module.

While `.apt` is more and more getting accepted, traditional text and binary file formats
are still commonly used in daily atom probe research practice. Not for all of these
formats formal specifications exist. This makes working with these formats
in software tools other than from AMETEK/Cameca trickier and error-prone.

A practical solution to raise at least awareness of this problem has been that scientists
collect examples (instances) of files in respective formats. Pieces of information about the
content and formatting of atom probe file formats were reported in the literature
(e.g. in the books by [D. Larson et al.](https://doi.org/10.1007/978-1-4614-8721-0) or
[B. Gault et al.](https://dx.doi.org/10.1007/978-1-4614-3436-8). Atom probers like Daniel Haley
have contributed substantially through raising awareness of the issue within the community.
Consequently, individuals of the community invested into reverse engineering efforts about
what these formats store and how this can be parsed using open-source software that is
developed within the atom probe community and beyond.

The `ifes_apt_tc_data_modeling` library bundles this knowledge highlighting though also
that there are still gaps in our understanding. From an academic point of view
these should be closed so that whenever possible atom probe data and metadata can
be always communicated clearly with respect to what do certain numbers mean, i.e.,
what are the semantics and concepts behind the numbers and data items. 

As an example, the `.pos` file format stores a table of number quadruples which mostly
are interpreted as reconstructed position and mass-to-charge-state ratio values.
Often the latter column is hijacked though to report conceptually different quantities
like identifier used to distinguish clusters of atoms. Which specific input data were
used, which parameterization was used for the reconstruction algorithm whose results
were stored in that `.pos` file. These questions pertaining to the workflow and
provenance along the data lifecycle remain unaddressed. Other technical issues exist
with file formats like `.pos` and `.epos`: These do not provide a [magic number](https://en.wikipedia.org/wiki/List_of_file_signatures)
that identifies the file as a true `.pos` file such that software tools and humans
could make substantiated assumptions.

Needs for improvement exist also for ranging definitions file formats like the commonly used `.rrng`, `.rng`,
and `.env` formats: These merely store the resulting ranging definitions but do not store details based on which
peak finding algorithm or even which mass-to-charge-state-ratio value array they were defined with.
[A more detailed discussion of these limitations is provided in the literature](https://doi.org/10.1017/S1431927621012241).

The `ifes_apt_tc_data_modeling` library was developed after observing that
many researchers in atom probe uses custom written code for reading atom probe
data via classical file formats. While for several formats this is a rather
simple programming exercise, it led though to parallel developments and many
implementations that target only specific use cases instead of a general enough
implementations with functionalities for all possible elements, ion types, and
edge cases.

