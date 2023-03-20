# atomprobe-data-modeling

## Mission:
Foster exchange about data models and work towards specifications
of file formats from the research field of atom probe microscopy.

## Documentation of file formats and data models in atom probe status quo
Detailed technical specifications of the file formats and data models are not
available for all formats in the field of atom probe microscopy.
A practical solution to address this limitation so far has been that scientists collect
examples in respective formats, so-called instances, and inspect and share these.
Based on this reverse engineering, individual atom probers have contributed to formulate
what can be considered likely candidates of the specification for several file formats.
This worked especially well for the POS and ePOS formats.

Pieces of information about file formats are reported in the literature in e.g.
books by D. Larson et al. and B. Gault et al. Atom probers like D. Haley have contributed
substantially to make the community aware of existent limitations and these
reverse engineering practices. AMETEK/Cameca is the key technology partner in atom probe.
The company has created an open file format called APT which improves the accessibility of
specific numerical data and some metadata. Individuals like M. Kühbach have driven the
implementation and communication of parsers for this APT file format.

Nowadays there is an increased interest and demand placed on atom probers by the funding
agencies that researchers should and have to make their research data management and data
stewardship better matching and more completely aligned to the aims and practices of
F.A.I.R. research principles. Therefore, it is useful to exchange more details about
data models and file formats as otherwise it is not foreseeable how atom probe data
can be made really interoperable with electronic lab notebooks, research data management
systems (RDMS), and related software tools.

In light of these challenges, the idea of understanding formats just by examples,
showed to be a slow and error-prone route as e.g. source code and workflows which have
been useed to write such files, and the associated input, workflow, and provenance information
might have been or not captured or/and the specific software tool(s) used might not have been
shared or made accessible for reviewing and analyzing their functions.

## Benefit and Next Steps
You can easily imagine that the more people will support us with this work the
more complete our understanding and knowledge about the available file formats
in atom probe microscopy can and will become. This can help all of us in the
long run to build software tools which are more reliable and thrustworthy and
also technically more robust when it comes to parse research data - irrespective
from which tools these come or how you would like to used these data further
throughout your journey of digitalizing your atom probe research.

The Python parsers in this repository are meant as a motivation to offer
immediate benefit to users. The collection of examples and technical discussions
via issues serves the more long-term aim. This is to arrive at a detailed technical
specification rather than more robust parsers only so that atom probe data can be
exchanged across tools irrespective of their formatting.

# Getting started as developer

You should create a virtual environment. We tested on Ubuntu with Python 3.7.
Newer versions of Python should work similarly (for this change 3.7) to the
desired version tag.

If you don't have Python 3.7 installed on your computer, follow these commands:
```
sudo add-apt-repository ppa:deadsnakes/ppa
sudo apt install python3.7 python3-dev libpython3.7-dev python3.7-venv
```

In some cases when using Python3.7 it was necessary to install python-numpy in others
not. So consider this if you run into issues when continuing this manual.
The following steps will install the ifes_apt_tc_data_modeling module in the
latest version.

```
mkdir <your-brand-new-folder>
cd <your-brand-new-folder>
pip install virtualenv
virtualenv --python=python3.7 .py37
source .py37/bin/activate

git clone git@github.com:atomprobe-tc/ifes_apt_tc_data_modeling.git
cd ifes_apt_tc_data_modeling
pip install .
pip install --upgrade pip
pip list
```

You can find instructions about how to use this tool in the tests/data jupyter notebook.
This notebook can be started from the command line inside the ifes_apt_tc_data_modeling
directory simply by calling.

```
jupyter-lab
```


## Support us with this work
Thank you very much for supporting this activity and your time.

## Feedback, questions
Feel free to drop us a message via creating an issue or commenting on one.
Feel invited to use the resources in this repository.

## Where to place your examples?
There is a *examples_with_provenance* and *examples_without_provenance*
sub-directory for each file format.

When you do know with which software and measured dataset you have created a file,
you should share the file and these pieces of information (software version). Do so by
naming at least the respective raw files. Ideally you share the examples via offering
a link to an external data repository such as Zenodo or other providers. This not only
avoids that this repository would get too much filled up with binary data.
Also it enables you to share clearly under which license you would like make your
example(s) accessible.

## Provenance if possible, plain examples if in doubt
Use the *examples_with_provenance* sub-directory. With this it is at least possible
to reproduce the file creation. A practical solution is to share (by uploading)
the screenshot of the complete IVAS/APSuite version info screen, including
the APSuite version, the CernRoot version, the CamecaRoot version, and the versions
of libraries used by APSuite. This can help other atom probers and AMETEK/Cameca
to improve their software as it enables them to identify inconsistencies
or bugs eventually easier.

Atom probers should be aware that file formats like POS, ePOS, or APT are neither
raw data nor follow a clear technical documentation. Therefore, all current file
formats are not meeting the FAIR principles.ey specified. Instead, refer to RRAW,
STR, RHIT and/or HITS files. Ideally, you add unique identifiers (such as SHA256
checksums) for each raw file. A documentation how you can do this was issued by 
your IFES APT TC colleagues [(How to hash your data)](https://github.com/oxfordAPT/hashlist).

If you cannot provide such detailed pieces of information to your work you can
still participate and support us a lot if you share your knowledge by adding at
least a link to a repository or file share with content in the relevant atom-probe
-specific file formats.
In this case, please use the *examples_without_provenance* directory.
While these examples are stripped of the context in which they were created
and used (provenance information), these examples can still be very useful
to run the file formats parsers against to make the parsers more robust.

# Background information
File formats, data models, in (almost every) research field may not be 
fully documented. A checklist of the necessary pieces of information and 
documentation required to call a data model, data schema, and/or file format
fully documented in accordance with the FAIR data and research software stewardship
principles is given below:

1. Each piece of information (bit/byte) is documented.
2. This documentation fulfills the FAIR principles, i.e.
   [Wilkinson et al., 2016](https://doi.org/10.1038/sdata.2016.18) and
   [Barker et al., 2022](https://doi.org/10.1038/s41597-022-01710-x)
   For binary files, tools like [kaitai struct](https://kaitai.io/) offer a
   solution to describe the exact binary information content in a data
   item. This can be a file but also the storage of a database entry or the
   response of a call to an API.  Let alone the binary structure is insufficient,
   tough.
3. To each piece of information there has to exist also a parameterized description,
   what this piece of information conceptually means. One way to arrive at such
   description is to use a data schema or ontology.
   It is important to mention that the concepts in this schema/ontology have
   unique identifier so that each data item/piece of information is identifiable
   as an instance of an entry in a database or a knowledge graph.
   This holds independently of which research data management system
   or electronic lab notebook is used.
4. In addition, it is very useful that timestamps are associated with 
   each data item (ISO8061 with time zone information) so that it is possible
   to create a timeline of the context in which and when the e.g. file was created.

The first and second point is known as a specification, while the third and fourth
point emphasize that the contextualization and provenance is key to make a
specification complete and useful.
