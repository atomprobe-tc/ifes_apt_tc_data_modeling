# Provenance

## Provenance and how to improve on it

Atom probers should be aware that file formats like `.pos`, `.epos`, or `.apt` are neither
raw data nor follow a clear technical documentation that is completely available to the public.
Therefore, all current file formats are not meeting the FAIR principles.

Instead, share `.rraw`, `.str`, `.rhit`, and `.hits` files when working with AMETEK/Cameca instruments.
Ideally, you add unique identifiers (such as SHA256 checksums) for each file.
A documentation how you can do this was issued by your IFES APT TC colleagues
[(How to hash your data)](https://github.com/oxfordAPT/hashlist).

