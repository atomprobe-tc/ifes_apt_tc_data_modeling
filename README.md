# atomprobe-data-modeling

## Mission:
Foster exchange about data models and work towards specifications
of file formats from the research field of atom probe microscopy.

# Getting started
You should create a virtual environment. We tested on Ubuntu with Python 3.8 and newer version.
In what follows the version (tag) 3.8 is a placeholder whereby we show how to proceed when using
Python 3.8. Using newer versions of Python should work the same by replacing 3.8 with the respective
version (tag).

Older versions of Python like 3.8 and 3.9 are available e.g. via the deadsnakes repository or via
conda. For using deadsnakes proceed with the following commands:
```
sudo add-apt-repository ppa:deadsnakes/ppa
sudo apt install python3.8 python3-dev libpython3.8-dev python3.8-venv
```

In some cases when using Python3.8, it was necessary to install python-numpy.
Please consider this if you run into issues when continuing with this manual.
The following steps will install the ifes_apt_tc_data_modeling module in the
latest version.

```
mkdir <your-brand-new-folder>
cd <your-brand-new-folder>
pip install virtualenv
virtualenv --python=python3.8 .py38
source .py38/bin/activate

git clone git@github.com:atomprobe-tc/ifes_apt_tc_data_modeling.git
cd ifes_apt_tc_data_modeling
python -m pip install --upgrade pip
python -m pip install pip-tools
python -m pip install -e .
python -m pip list
```

## Additional steps to do when working with jupyter notebooks
By default the functionalities are offered as a library for Python programmers.
For developers and users who would like to try using the library a convenient
way is via jupyter notebooks. You can find instructions about how to use this tool
in the tests/data jupyter notebook. This notebook can be started from the command
line inside the ifes_apt_tc_data_modeling directory simply by calling.
If you would like to use a jupyter notebook jupyter has to be installed as
it will not be installed by default. To achieve this perform the following actions:

```
python -m pip install -e ".[dev]"
python -m pip list

jupyter-lab
```

## Documentation of file formats and data models in atom probe status quo
Detailed technical specifications of the file formats and data models are not available for
most formats in the field of atom probe microscopy. A practical solution to address this 
limitation has been so far that scientists collect example files formatted in respective formats.

These so-called instances were inspected and shared with colleagues. In summary, individual
atom probers have contributed to formulate what can be considered likely candidates
of specifications for several file formats via reverse engineering.
This worked especially well for the POS and ePOS formats.

Pieces of information about file formats were reported in the literature (e.g.
the books by D. Larson et al. and B. Gault et al.). Atom probers like D. Haley have contributed
substantially to make the community aware of existent limitations and these reverse engineering
practices. AMETEK/Cameca is the key technology partner in atom probe. They have developed
an open file format called APT which improves the accessibility of specific numerical data and
some metadata. Individuals like M. Kühbach have driven the implementation and communication of
parsers for this APT file format.

Nowadays there is an increased interest and demand placed on atom probers by the funding agencies
that researchers should or even have to make their research data management and data stewardship
better matching and more completely aligned to the aims and practices of the F.A.I.R.
principles of data stewardship. Therefore, it is useful to exchange more details about
data models and file formats. Otherwise, it is not foreseeable how atom probe data can be made
really interoperable with electronic lab notebooks, research data management
systems (RDMS), and related software tools for data analyses.

In light of these challenges, the idea of understanding formats just by examples, showed to be a
slow and error-prone route as e.g. source code and workflows which have been useed to write such
files, and the associated input, workflow, and provenance information has typically not been captured.
Or the specific software tool(s) used might not have been shared or made accessible for review
by the atom probe community.

## Benefit and Next Steps
You can easily imagine that the more people will support this work the more complete a public
understanding and knowledge about the available file formats in atom probe microscopy will become.
This can help all of us in the long run to build software tools which are more reliable, yield
thrustworthy results and are technically more robust when it comes to parsing research data.
Irrespective from which tools these data and metadata come or how one would like to used these data.

The Python parsers in this repository are meant as a motivation to offer immediate benefit for users.
The collection of examples and technical discussions via issues serves the more long-term aim.
This is to arrive at a detailed technical specification rather than having more robust parsers only
so that atom probe data can be exchanged across tools irrespective of their formatting.

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

# Background information
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
