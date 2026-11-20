# chargefw-ccd

The atoms and bonds of every component of the wwPDB
[Chemical Component Dictionary](https://www.wwpdb.org/data/ccd), packaged as one data file for
[ChargeFW](https://chargefw.readthedocs.io/).

```
pip install chargefw-ccd
```

```python
import chargefw_ccd

chargefw_ccd.path()  # location of the data file
chargefw_ccd.__version__  # for example 2026.11.19
```

The version is the date of the most recently modified component in the dictionary the package was
built from.

The data come from `components.cif` of the wwPDB, which is available under
[CC0 1.0](https://creativecommons.org/publicdomain/zero/1.0/).
