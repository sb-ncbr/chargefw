# Write charges to RDKit molecules

This script reads an SDF file with RDKit, calculates charges for all records, and writes them back as
RDKit atom properties. For [paracetamol.sdf](../getting-started/paracetamol.sdf):

```console
$ pip install --pre "chargefw[rdkit]"
$ python charge_rdkit_sdf.py paracetamol.sdf charged.sdf --method eem
Model: eem/EEM_NEEMP_ccd2016_npa
'paracetamol': total=+0.000 most positive=C2 (+0.478)
```

The output SDF stores the charges in an `atom.dprop.ChargeFWPartialCharge` property, which RDKit turns
back into the `ChargeFWPartialCharge` property of each atom when it reads the file. The method and
parameter set are stored as the molecule properties `ChargeFWMethod` and `ChargeFWParameterSet`:

```text
>  <atom.dprop.ChargeFWPartialCharge>  (1)
-0.59586795804550241 0.47829714521921546 -0.63115993457608788 ...

>  <ChargeFWMethod>  (1)
eem
```

Records that RDKit cannot read are reported and skipped. The other records must all have explicit
hydrogens. The script kekulizes aromatic rings so that ChargeFW gets single and double bonds. All
records are calculated with one model, as described in
[One model for all molecules](../interfaces/python.md#one-model-for-all-molecules); pass
`--parameter-set` to choose it.

Download {download}`charge_rdkit_sdf.py`.

```{literalinclude} charge_rdkit_sdf.py
:language: python
```
