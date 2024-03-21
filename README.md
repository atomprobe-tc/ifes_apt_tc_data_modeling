# ifes_apt_tc_data_modeling

## Mission:
Foster exchange in the community of atom probe research to exchange about and
document information content and formatting in their research field.
Work towards ideally semantically specified file formats and data models.

## Getting started:

### Create an environment
To use this library create a conda or a virtual environment. We tested on Ubuntu with Python 3.8 and newer version.
In what follows the version (tag) 3.8 is a placeholder whereby we show how to proceed when using e.g. Python version 3.8.
Using newer versions of Python should work the same by replacing 3.8 with the respective version (tag). As of 2024,
using Python in versions higher than 3.9 becomes more and more common. The support for users to install modern
Python version has also improved. Therefore, the following commands typically enable you to create a 
specifically-versioned virtual environment:

```
mkdir <your-brand-new-folder>
cd <your-brand-new-folder>
pip install virtualenv
virtualenv --python=python3.8 .py38
source .py38/bin/activate
```

If you wish to use or still demand to use older versions of Python, like 3.8 or 3.9, you can conveniently install them
via the deadsnakes repository (or via conda). For using deadsnakes proceed with the following commands:

```
sudo add-apt-repository ppa:deadsnakes/ppa
sudo apt install python3.8 python3-dev libpython3.8-dev python3.8-venv
```

In some cases when using Python3.8, it was necessary to install python-numpy.
Please consider this if you run into issues when continuing with this manual.

### Install the ifes_apt_tc_data_modeling modules as a user

```
git clone git@github.com:atomprobe-tc/ifes_apt_tc_data_modeling.git
cd ifes_apt_tc_data_modeling
python -m pip install --upgrade pip
python -m pip install pip-tools
python -m pip install -e .
python -m pip install -e ".[dev]"
python -m pip list
```

### Additional steps to perform when you are a developer

```
python -m pip install -e ".[dev]"
python -m pip list
jupyter-lab
```

## Context, status quo, file formats used for atom probe research
A lack of detailed technical specifications of the file formats and a lack of usage of magic numbers as identifiers for specific file formats
are a key blocker to parsing and semantic interpretation of information content stored in current file formats within the research field of
atom probe microscopy.

A practical solution to raise at least awareness of this problem has been that scientists collect examples (instances)
of files in respective formats. Pieces of information about the content and formatting of atom probe file formats were reported in the literature
(e.g. in the books by D. Larson et al. https://doi.org/10.1007/978-1-4614-8721-0 or B. Gault et al. http://dx.doi.org/10.1007/978-1-4614-3436-8 ).
Atom probers like D. Haley have contributed substantially through raising awareness of the issue within the community.

AMETEK/Cameca is the key technology partner in atom probe. AMETEK has developed an open file format called APT which has improved
the accessibility of specific numerical data and some metadata. Individuals like M. Kühbach have driven the implementation and
communication of parsers for this APT file format. There are ongoing efforts by both AMETEK and the scientific community to extent the APT file format
with additional metadata. The main motivation behind these newer efforts is to improve the interoperability between research data collected
within the IVAS/APSuite software and third-party software including research data management systems.
Currently, most metadata have to be entered manually via e.g. electronic lab notebooks if one were to use or register atom probe
data in solutions other than those developed by AMETEK.

Nowadays, there is a global desire, a push by research funding agencies, and an increased interest of atom probers
to make their research data and knowledge generation process better matching and more completely aligned to the aims
and practices of the F.A.I.R. principles of research data stewardship and FAIR4RS research software development.

Therefore, it is useful to exchange more details about data models and file formats. Otherwise, it is not foreseeable
how atom probe data can be made really interoperable with electronic lab notebooks, research data management
systems (RDMS), and related software tools for data analyses, especially not if these tools should ever work with
solutions from the stack of semantic web technologies. We are convinced there are substantial opportunities with making
atom probe research communication more substantiated, the research itself better reproducible, and with enabling
automated contextualization of atom probe research via computational agents.

In light of these challenges, the idea of understanding formats just by examples, showed to be a slow and error-prone route
as e.g. source code and workflows which have been used to write such files lack provenance information. As an example,
the POS files only store a table of number quadruples which mostly are interpreted as reconstructed position and mass-to-charge-
state ratio values but often are hijacked to report conceptually different quantities like identifier used to distinguish clusters of
atoms. Nowhere in a POS file a magic number could identify the file as to be truely a POS file and no something else based on
which software tools and human could make a substantiated assumption. Nowhere does the POS file document from which
content and which tools it was generated. The situation is currently still similarly poor for ranging definitions files such as RRNG, RNG,
or ENV: These merely store the resulting ranging definitions but no details based on which peak finding algorithm or even which
mass-to-charge-state-ratio value array they were defined with. M. Kühbach et al. have summarized a more detailed discussion
about these limitations https://doi.org/10.1017/S1431927621012241.

## How can you support this work?
As a user with contacting us and providing examples of file formats. As a member of a company by documenting your file format
and getting in contact to work together on improving the situation. Thank you very much for supporting this activity and your time.

## Feedback, questions
Feel free to drop us a message via creating an issue or commenting on one. 

## Background information
File formats, data models, in (almost every) research field may not be fully documented.
A checklist of the necessary pieces of information and documentation required to call a
data model, data schema, and/or file format fully documented in accordance with the
FAIR data and research software stewardship principles is given below:

1. Each piece of information (bit/byte) is documented.
2. This documentation fulfills the FAIR principles, i.e.
   [Wilkinson et al., 2016](https://doi.org/10.1038/sdata.2016.18) and
   [Barker et al., 2022](https://doi.org/10.1038/s41597-022-01710-x)
   For binary files, tools like [kaitai struct](https://kaitai.io/) offer a
   solution to describe the exact binary information content in a data
   item. This can be a file but also the storage of a database entry or the
   response of a call to an API.
   Let alone the binary structure is insufficient tough.
3. To each piece of information there has to exist also a parameterized description,
   what this piece of information conceptually means. One way to arrive at such
   description is to use a data schema or ontology.
   It is important to mention that the concepts in this schema/ontology have
   unique identifier so that each data item/piece of information is identifiable
   as an instance of an entry in a database or a knowledge graph.
   This holds independently of which research data management system
   or electronic lab notebook is used.
4. In addition, it is very useful if timestamps are associated with each data item
   (ISO8061 including time zone information) so that it is possible to create a
   timeline of the context in which and when the e.g. file was created.

The first and second point is known as a specification, while the third and fourth
point emphasize that the contextualization and provenance is key to make a
specification complete and useful.

## Where to place your examples?
There is a *examples_with_provenance* and *examples_without_provenance* sub-directory for each file format.

When you do know with which software and measured dataset you have created a file,
you should share the file and these pieces of information (software version). Do so by
naming at least the respective raw files. Ideally, you share the examples via offering
a link to an external data repository such as Zenodo or other providers. This not only
avoids that this repository would get too much filled up with binary data.
Also it enables you to share clearly under which license you would like make your
example(s) accessible.

## Provenance if possible, plain examples if in doubt
Use the *examples_with_provenance* sub-directory. With this it is at least possible
to reproduce the file creation. A practical solution is to share (by uploading)
the screenshot of the complete IVAS/APSuite version info screen, including
the APSuite version, the CERN Root version, the CamecaRoot version, and the versions
of libraries used by APSuite. This can help other atom probers and AMETEK/Cameca
to improve their software as it will enable them to identify inconsistencies.

Atom probers should be aware that file formats like POS, ePOS, or APT are neither
raw data nor follow a clear technical documentation. Therefore, all current file
formats are not meeting the FAIR principles. Instead, share RRAW, STR, RHIT, and HITS files.
Ideally, you add unique identifiers (such as SHA256 checksums) for each file.
A documentation how you can do this was issued by your IFES APT TC colleagues
[(How to hash your data)](https://github.com/oxfordAPT/hashlist).

If you cannot provide such detailed pieces of information, you can still participate
and support us a lot if you share your knowledge by adding at least a link to a repository
or file share with content in the relevant atom-probe-specific file formats.

In this case, please use the *examples_without_provenance* directory.
While these examples are stripped of the context in which they were created
and used (provenance information), these examples can still be very useful
to run the file formats parsers against to make the parsers more robust, i.e.
that these can pick up formatting issues and act accordingly.
