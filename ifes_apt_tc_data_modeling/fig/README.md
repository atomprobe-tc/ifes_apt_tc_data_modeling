## Getting started
Matlab figures are not necessarily directly parsable using Python tools.
Therefore, atom probers who store their ranging definitions as Matlab
figures should first run the following Matlab script `matlab/matlab_fig_to_txt.m`
to extract the ranging definitions into a text file.

The text file is very simple. Each ranging definitions is a single line
with three parts. The first part is a human-readable description of the
ion (element, isotope, molecular ion). We follow the naming convention
of P. Felfer's atom-probe-toolbox. Examples are shown below.
The second part is the left (minimum) while the third
part is the right (maximum bound) bound of the mass-to-charge-state
ratio interval identified as the ion named in the first part.

The file must not contain other definitions

```
16O 1H2+ 1.796520e+01 1.808020e+01
12C2+ 2.394526e+01 2.407026e+01
50Cr++ 2.494527e+01 2.503027e+01
```
