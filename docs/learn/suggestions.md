# Specifications

## Improve on building and communicating file format specifications

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




