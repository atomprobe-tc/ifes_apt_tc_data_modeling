## An example for reading XML-serialized program state with data content exported from IVAS prior AP Suite

Many atom probe groups have still the resultant XML-serialized state files within their projects which IVAS
prior AP Suite - to the best of my knowledge - generated to document data and metadata within the GUI.

This imago reader is a very simple example for a function showing how one can extract information from these
XML files. The example substantiates clearly why modern approaches for storing metadata data like
the NeXus proposal for atom probe (www.github.com/FAIRmat-NFDI/nexus_definitions) are much cleaner
and likely easier for end users to interact with.

Indeed, the here exemplified serialized XML file is in a format that when represented in python mixes
recurrent nested lists within dictionaries, a signature of the fact how Java-based applications like older
IVAS versions serialize information. Such a data structure is not directly flattenable and hence many
checks are required to fish even the simplest pieces of information.

I expect that this parser does not work out of the box for many examples given that the formatting of
object-serialized content depends heavily on the specific implementation of the host application (IVAS)
and thus the IVAS version. The main idea behind the example is to show that information can be extracted
though which could be useful for legacy purposes. Also I would like to motivate that nowadays one should
use data structures that are more conveniently parsable like aforementioned NeXus data schema and
respective HDF5 implementations as e.g. exemplified in the www.gitlab.com/paraprobe/paraprobe-toolbox
reference implementation for the atom probe community.
