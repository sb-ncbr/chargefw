# Methods and parameter sets

A method is an empirical charge model, such as EEM or PEOE. Most methods also need a parameter set: the
element-specific values that the authors of a method, or later studies, fitted to reproduce a chosen
kind of reference charges, for example B3LYP/6-31G* NPA charges. ChargeFW includes 20 methods; the
[parameter-set reference](../reference/parameter-sets.md) lists the sets with
their publications. All methods calculate charges of finite molecules; periodic systems such as
crystals are outside their scope.

## Choosing a method

If you know which reference charges you need, choose the method and the parameter set fitted to them,
and name both. On the command line:

```bash
chargefw calculate --method eem --parameter-set EEM_Cheminf_b3lyp_npa molecule.sdf output
```

In Python:

```python
result = chargefw.calculate(molecules, method="eem", parameter_set="EEM_Cheminf_b3lyp_npa")
```

ChargeFW calculates exactly what you name. If the method or the parameter set cannot handle your
molecules, the calculation fails and reports why; ChargeFW does not switch to another model.

If you name nothing, or only the method, ChargeFW chooses for you. It tries the methods in the order of
the [table below](#available-methods) and, for each method, its parameter sets in the order of the
parameter-set reference. It uses the first combination that applies to every molecule in the
calculation, so all molecules get charges from the same model. The result records which method and
parameter set were used.

To see why a method or parameter set was not used, without calculating anything:

```bash
chargefw applicability molecule.sdf
```

In Python, `chargefw.assess()` returns the runnable plans and, in `rejections`, the reasons for every
rejected candidate; see [Python](../interfaces/python.md#advanced-assessment-and-plan-reuse).

A new ChargeFW release can add parameter sets and so change what automatic selection picks. Name the
method and the parameter set whenever results must stay the same between releases.

## Available methods

The table lists the methods in the order automatic selection tries them. **Coordinates** shows whether
the method needs 3D coordinates; methods that do produce one set of charges per conformer.
**Parameter sets** gives the number of bundled sets; methods marked "built-in" use element properties
included in ChargeFW instead. **Cutoff/cover** marks methods that can run the approximations for
[large structures](large-structures.md), and **Fixed ions** those that accept
[fixed ions](fixed-ions.md). **Time** is the calculation cost for $n$ atoms and $m$ bonds.

| Method | ID | Coordinates | Parameter sets | Cutoff/cover | Fixed ions | Time |
| --- | --- | --- | --- | --- | --- | --- |
| [SQE+qp](#sqeqp) | `sqeqp` | yes | 3 | yes | yes | $O(m^3)$ |
| [EEM](#eem) | `eem` | yes | 28 | yes | yes | $O(n^3)$ |
| [ABEEM](#abeem) | `abeem` | yes | 1 | yes | | $O((n+m)^3)$ |
| [SFKEEM](#sfkeem) | `sfkeem` | yes | 1 | yes | | $O(n^3)$ |
| [QEq](#qeq) | `qeq` | yes | 1 | yes | | $O(n^3)$ |
| [SMP/QEq](#smpqeq) | `smpqeq` | yes | 1 | | | $O(n^3)$ |
| [EQeq](#eqeq) | `eqeq` | yes | built-in | yes | | $O(n^3)$ |
| [EQeq+C](#eqeqc) | `eqeqc` | yes | 2 | yes | | $O(n^3)$ |
| [DelRe](#delre) | `delre` | | 1 | | | $O(n^3)$ |
| [PEOE](#peoe) | `peoe` | | 1 | | | $O(n + m)$ |
| [MPEOE](#mpeoe) | `mpeoe` | | 1 | | | $O(n + m)$ |
| [GDAC](#gdac) | `gdac` | yes | 1 | | | $O(n + m)$ |
| [SQE](#sqe) | `sqe` | yes | 3 | yes | yes | $O(m^3)$ |
| [SQE+q0](#sqeq0) | `sqeq0` | yes | 3 | yes | yes | $O(m^3)$ |
| [MGC](#mgc) | `mgc` | | built-in | | | $O(n^3)$ |
| [KCM](#kcm) | `kcm` | | 2 | | | $O(n^3)$ |
| [TSEF](#tsef) | `tsef` | | 1 | | | $O(n^3)$ |
| [DENR](#denr) | `denr` | | 1 | | | $O(n^3)$ |
| [Charge2](#charge2) | `charge2` | | 1 | | | $O(n + m)$ |
| [VEEM](#veem) | `veem` | | built-in | | | $O(n)$ |

Two further IDs are not charge models: `formal` returns the formal charges of your input unchanged, and
`dummy` returns zero for every atom. Automatic selection tries them last, so they are used only when no
charge model applies.

`chargefw methods METHOD` on the command line, or `chargefw.methods["METHOD"]` in Python, shows a
method's options and their defaults.

## Method descriptions

The symbols used throughout: $q_i$ stands for the charge on the $i$th atom, $Q$ is the total charge of
the molecule, $N$ is the number of atoms, and $R_{i,j}$ is the distance between atoms $i$ and $j$ in Å.

Methods derived from one another are kept together. Such families follow the year of their first
method, and the methods within a family follow their own years.

The figure shows the families by year of publication; a line joins a method to the one it derives
from. Select a method to go to its description.

<svg viewBox="0 0 800 543" role="img" aria-label="Families of the methods by year of publication" style="width: 100%; max-width: 800px; font-family: inherit">
<line x1="60" y1="498" x2="790" y2="498" style="stroke: var(--color-background-border, #ccc)"/>
<text x="10" y="503" style="fill: var(--color-foreground-muted, #666); font-size: 15px">1960</text>
<line x1="60" y1="428" x2="790" y2="428" style="stroke: var(--color-background-border, #ccc)"/>
<text x="10" y="433" style="fill: var(--color-foreground-muted, #666); font-size: 15px">1970</text>
<line x1="60" y1="358" x2="790" y2="358" style="stroke: var(--color-background-border, #ccc)"/>
<text x="10" y="363" style="fill: var(--color-foreground-muted, #666); font-size: 15px">1980</text>
<line x1="60" y1="288" x2="790" y2="288" style="stroke: var(--color-background-border, #ccc)"/>
<text x="10" y="293" style="fill: var(--color-foreground-muted, #666); font-size: 15px">1990</text>
<line x1="60" y1="218" x2="790" y2="218" style="stroke: var(--color-background-border, #ccc)"/>
<text x="10" y="223" style="fill: var(--color-foreground-muted, #666); font-size: 15px">2000</text>
<line x1="60" y1="148" x2="790" y2="148" style="stroke: var(--color-background-border, #ccc)"/>
<text x="10" y="153" style="fill: var(--color-foreground-muted, #666); font-size: 15px">2010</text>
<line x1="60" y1="78" x2="790" y2="78" style="stroke: var(--color-background-border, #ccc)"/>
<text x="10" y="83" style="fill: var(--color-foreground-muted, #666); font-size: 15px">2020</text>
<line x1="355" y1="20" x2="355" y2="538" style="stroke: var(--color-foreground-muted, #666); stroke-dasharray: 4 4"/>
<text x="205" y="22" text-anchor="middle" style="fill: var(--color-foreground-primary, #000); font-size: 16px; font-weight: bold">Topology</text>
<text x="580" y="22" text-anchor="middle" style="fill: var(--color-foreground-primary, #000); font-size: 16px; font-weight: bold">Geometry</text>
<line x1="250" y1="294" x2="200" y2="343" style="stroke: var(--color-foreground-muted, #666); stroke-width: 1.2"/>
<line x1="420" y1="217" x2="250" y2="273" style="stroke: var(--color-foreground-muted, #666); stroke-width: 1.2"/>
<line x1="100" y1="168" x2="110" y2="196" style="stroke: var(--color-foreground-muted, #666); stroke-width: 1.2"/>
<line x1="470" y1="245" x2="520" y2="301" style="stroke: var(--color-foreground-muted, #666); stroke-width: 1.2"/>
<line x1="560" y1="182" x2="520" y2="301" style="stroke: var(--color-foreground-muted, #666); stroke-width: 1.2"/>
<line x1="600" y1="161" x2="680" y2="266" style="stroke: var(--color-foreground-muted, #666); stroke-width: 1.2"/>
<line x1="690" y1="140" x2="680" y2="266" style="stroke: var(--color-foreground-muted, #666); stroke-width: 1.2"/>
<line x1="642" y1="108" x2="668" y2="129" style="stroke: var(--color-foreground-muted, #666); stroke-width: 1.2"/>
<line x1="755" y1="182" x2="680" y2="266" style="stroke: var(--color-foreground-muted, #666); stroke-width: 1.2"/>
<line x1="755" y1="140" x2="755" y2="161" style="stroke: var(--color-foreground-muted, #666); stroke-width: 1.2"/>
<line x1="755" y1="77" x2="755" y2="119" style="stroke: var(--color-foreground-muted, #666); stroke-width: 1.2"/>
<a href="#delre"><text x="150" y="513" text-anchor="middle" style="fill: var(--color-foreground-primary, #000); font-size: 15px; font-weight: bold">DelRe<title>DelRe, 1958</title></text></a>
<a href="#peoe"><text x="200" y="359" text-anchor="middle" style="fill: var(--color-foreground-primary, #000); font-size: 15px; font-weight: bold">PEOE<title>PEOE, 1980</title></text></a>
<a href="#charge2"><text x="100" y="345" text-anchor="middle" style="fill: var(--color-foreground-primary, #000); font-size: 15px; font-weight: bold">Charge2<title>Charge2, 1982</title></text></a>
<a href="#mpeoe"><text x="250" y="289" text-anchor="middle" style="fill: var(--color-foreground-primary, #000); font-size: 15px; font-weight: bold">MPEOE<title>MPEOE, 1990</title></text></a>
<a href="#mgc"><text x="110" y="212" text-anchor="middle" style="fill: var(--color-foreground-primary, #000); font-size: 15px; font-weight: bold">MGC<title>MGC, 2001</title></text></a>
<a href="#kcm"><text x="100" y="163" text-anchor="middle" style="fill: var(--color-foreground-primary, #000); font-size: 15px; font-weight: bold">KCM<title>KCM, 2008</title></text></a>
<a href="#tsef"><text x="190" y="163" text-anchor="middle" style="fill: var(--color-foreground-primary, #000); font-size: 15px; font-weight: bold">TSEF<title>TSEF, 2008</title></text></a>
<a href="#denr"><text x="280" y="163" text-anchor="middle" style="fill: var(--color-foreground-primary, #000); font-size: 15px; font-weight: bold">DENR<title>DENR, 2008</title></text></a>
<a href="#veem"><text x="150" y="142" text-anchor="middle" style="fill: var(--color-foreground-primary, #000); font-size: 15px; font-weight: bold">VEEM<title>VEEM, 2011</title></text></a>
<a href="#gdac"><text x="420" y="212" text-anchor="middle" style="fill: var(--color-foreground-primary, #000); font-size: 15px; font-weight: bold">GDAC<title>GDAC, 2001</title></text></a>
<a href="#eem"><text x="520" y="317" text-anchor="middle" style="fill: var(--color-foreground-primary, #000); font-size: 15px; font-weight: bold">EEM<title>EEM, 1986</title></text></a>
<a href="#abeem"><text x="470" y="240" text-anchor="middle" style="fill: var(--color-foreground-primary, #000); font-size: 15px; font-weight: bold">ABEEM<title>ABEEM, 1997</title></text></a>
<a href="#sfkeem"><text x="560" y="177" text-anchor="middle" style="fill: var(--color-foreground-primary, #000); font-size: 15px; font-weight: bold">SFKEEM<title>SFKEEM, 2006</title></text></a>
<a href="#qeq"><text x="680" y="282" text-anchor="middle" style="fill: var(--color-foreground-primary, #000); font-size: 15px; font-weight: bold">QEq<title>QEq, 1991</title></text></a>
<a href="#smpqeq"><text x="600" y="156" text-anchor="middle" style="fill: var(--color-foreground-primary, #000); font-size: 15px; font-weight: bold">SMP/QEq<title>SMP/QEq, 2009</title></text></a>
<a href="#eqeq"><text x="690" y="135" text-anchor="middle" style="fill: var(--color-foreground-primary, #000); font-size: 15px; font-weight: bold">EQeq<title>EQeq, 2012</title></text></a>
<a href="#eqeqc"><text x="610" y="114" text-anchor="middle" style="fill: var(--color-foreground-primary, #000); font-size: 15px; font-weight: bold">EQeq+C<title>EQeq+C, 2015</title></text></a>
<a href="#sqe"><text x="755" y="177" text-anchor="middle" style="fill: var(--color-foreground-primary, #000); font-size: 15px; font-weight: bold">SQE<title>SQE, 2006</title></text></a>
<a href="#sqeq0"><text x="755" y="135" text-anchor="middle" style="fill: var(--color-foreground-primary, #000); font-size: 15px; font-weight: bold">SQE+q0<title>SQE+q0, 2012</title></text></a>
<a href="#sqeqp"><text x="755" y="72" text-anchor="middle" style="fill: var(--color-foreground-primary, #000); font-size: 15px; font-weight: bold">SQE+qp<title>SQE+qp, 2021</title></text></a>
</svg>

### Topology-based methods

These methods use only the molecular graph: atoms, bonds, and formal charges. GDAC, which derives
from PEOE and MPEOE, also needs coordinates.

(delre)=
#### DelRe

The method of Del Re ([Del Re 1958](https://doi.org/10.1039/JR9580004031)) starts with a system of linear
equations:

$$
\delta_i = \delta_i^0 + \sum_j \gamma_{i \leftarrow j}\,\delta_j
$$

where $j$ iterates over atoms bonded to $i$, $\delta_i^0$ is an atom parameter, and $\gamma$ is a bond
parameter with a separate value for each direction of the bond. Solving for $\delta$ gives the bond
charges:

$$
q_{i,j} = \frac{\delta_i - \delta_j}{2\varepsilon_{i,j}}
$$

where $\varepsilon_{i,j}$ is another bond parameter. Atom $i$ loses the bond charge and atom $j$ gains
it; the charge on each atom is the sum over its bonds.

Every connected component must be formally neutral.

(peoe)=
#### PEOE

Partial Equalization of Orbital Electronegativity
([Gasteiger 1980](https://doi.org/10.1016/0040-4020(80)80168-2)) is an iterative scheme in which charge
moves along the bonds from the more electropositive atom to the more electronegative one. The
electronegativity of an atom is a function of its charge and is recomputed at the start of each
iteration $\alpha$:

$$
\chi_i = A_i + B_i q_i + C_i q_i^2
$$

where $A$, $B$, and $C$ are atom parameters. For every bond $i$–$j$ with $\chi_j > \chi_i$, the charge
of atom $i$ increases and the charge of atom $j$ decreases by:

$$
\Delta q = \left(\frac{1}{2}\right)^{\alpha} \frac{\chi_j - \chi_i}{A_i + B_i + C_i}
$$

The damping factor $(1/2)^\alpha$ accounts for the electrostatic field of the charges already
transferred. For hydrogen, the denominator is a hydrogen-specific value. The final charges are the initial
charges plus the transfers of all iterations.

`iters` sets the number of iterations (default 6). With the default `initial_charges=zero`, every
connected component must be formally neutral; `initial_charges=formal` starts from the formal charges
and allows charged molecules.

(mpeoe)=
#### MPEOE

Modified Partial Equalization of Orbital Electronegativity
([No 1990](https://doi.org/10.1021/j100374a066)) differs from [PEOE](#peoe) by expressing electronegativity
as a linear function of charge:

$$
\chi_i = A_i + B_i q_i
$$

The damping factor $(1/2)^\alpha$ is replaced by a bond parameter $f_{i,j}$:

$$
\Delta q = f_{i,j}^{\,\alpha}\, \frac{\chi_j - \chi_i}{A_i + B_i}
$$

As in PEOE, hydrogen uses a hydrogen-specific denominator.

`iters` sets the number of iterations (default 7). `initial_charges` and the charge requirement are the
same as for PEOE.

(gdac)=
#### GDAC

Geometry-Dependent Atomic Charges ([Cho 2001](https://doi.org/10.1021/jp0023213)) modifies the damping
factor of [MPEOE](#mpeoe) to depend on the geometry:

$$
f_{i,j} = 1 - \frac{R_{i,j}}{R_i^\mathrm{vdW} + R_j^\mathrm{vdW}}
$$

where $R^\mathrm{vdW}$ is the van der Waals radius. The denominator is $A_i + B_i$ for all atoms,
including hydrogen.

`iters` sets the number of iterations (default 7). Every connected component must be formally neutral.

(charge2)=
#### Charge2

Charge2 ([Abraham 1982](https://doi.org/10.1002/jcc.540030316)) is an iterative method in which atoms one,
two, and three bonds apart contribute charge increments to a central atom. Each increment is transferred
between two atoms, so the total charge is conserved:

$$
\begin{aligned}
q_i(\alpha) &= \sum_j \frac{\chi_j - \chi_i}{a}\\
q_i(\beta) &= \sum_k \frac{(\chi_k - \chi_\mathrm{H}) P_i}{b}\\
q_i(\gamma) &= \sum_l \frac{(\chi_l - \chi_\mathrm{H}) P_i}{bc}
\end{aligned}
$$

where the sums run over all paths of one, two, and three bonds from atom $i$ that do not return to an
atom already on the path, and $j$, $k$, and $l$ are the atoms at their ends. In rings, an atom reached
by several such paths contributes once for each of them, and a bonded atom can also contribute through
a three-bond path. $\chi$ is an atom parameter, $\chi_\mathrm{H} = 7.17$ is the electronegativity of
hydrogen, and $b$ and $c$ are common parameters. The one-bond increment is taken from atom $j$; the two- and three-bond increments are taken
from the neighbor of $i$ on the path. The denominator $a$ is one of three common parameters, chosen by
whether both atoms are in period 2, one of them is hydrogen, or neither. The polarizability $P_i$
depends on the charge:

$$
P_i = \max\left(0,\; P_i^0 \left(1 + \alpha(q_i^0 - q_i)\right)\right)
$$

where $P^0$ and $q^0$ are atom parameters and $\alpha$ is a common parameter. The iteration stops when no
charge changes by more than $10^{-5}$.

`iters` sets the maximum number of iterations (default 10). Every connected component must be formally
neutral.

(mgc)=
#### MGC

Molecular Graph Charge model ([Oliferenko 2001](https://doi.org/10.1002/poc.378)) uses the molecular
graph, inspired by electrical circuits and Kirchhoff's laws; no atomic coordinates or parameters are
used. MGC constructs the matrix:

$$
\mathbf{S} = \mathbf{I} + \mathbf{D} - \mathbf{A}
$$

where $\mathbf{A}$ is the connectivity matrix with bond orders as entries, $\mathbf{D}$ is the diagonal
matrix of their row sums, and $\mathbf{I}$ is the identity matrix. Equalized electronegativities are
obtained from those of isolated atoms by solving:

$$
\mathbf{S} \boldsymbol{\chi} = \boldsymbol{\chi}^0
$$

The charges are the differences between equalized and isolated-atom electronegativities, divided by the
geometric mean $\chi_M$ of the isolated-atom electronegativities:

$$
\mathbf{q} = \frac{\boldsymbol{\chi} - \boldsymbol{\chi}^0}{\chi_M}
$$

Every connected component must be formally neutral.

(kcm)=
#### KCM

Kirchhoff Charge Model ([Yakovenko 2008](https://doi.org/10.1002/jcc.20892)) builds a Laplacian matrix
$\mathbf{L}$ as:

$$
\mathbf{L} = \mathbf{B}^{\mathsf T} \mathbf{W} \mathbf{B}
$$

where $\mathbf{B}$ is the bond–atom incidence matrix and $\mathbf{W}$ is a diagonal softness matrix with
$w_{i-j} = 1 / (\eta_i + \eta_j)$ for each bond, where $\eta$ is the atom hardness. The charges are:

$$
\mathbf{q} = (\mathbf{I} + \mathbf{L})^{-1} \boldsymbol{\chi}^0 - \boldsymbol{\chi}^0
$$

where $\boldsymbol{\chi}^0$ is the vector of atom electronegativities.

Every connected component must be formally neutral.

(denr)=
#### DENR

Dynamical Electronegativity Relaxation ([Shulga 2008](https://doi.org/10.1080/10629360701844142)) is an
iterative scheme in which charges are derived using a Laplacian matrix:

$$
\mathbf{q}^{(n+1)} = \left(\mathbf{I} + \Delta t\, \mathbf{L} \boldsymbol{\eta}\right)^{-1}
\left(\mathbf{q}^{(n)} - \Delta t\, \mathbf{L} \boldsymbol{\chi}^0\right)
$$

where $\mathbf{L}$ is the Laplacian of the molecular graph, $\boldsymbol{\eta}$ is the diagonal matrix
of atom hardnesses, $\boldsymbol{\chi}^0$ is the vector of atom electronegativities, and
$\mathbf{q}^{(0)} = \mathbf{0}$.

`step` sets $\Delta t$ (default 0.1) and `iters` the number of steps (default 3). Every connected
component must be formally neutral.

(tsef)=
#### TSEF

Topologically Symmetrical Energy Function ([Shulga 2008](https://doi.org/10.1080/10629360701844142)) is
based on electronegativity equalization, like [EEM](#eem), but replaces the Euclidean distance with the
bond distance, which makes TSEF independent of the conformation:

$$
\chi_i = \chi_i^0 + \eta_i q_i + \sum_{j \neq i} \phi_{i,j} q_j
$$

$$
\phi_{i,j} = 14.4 \cdot K(\mathrm{MDP}_{i,j}) \cdot \frac{1}{0.84\, \mathrm{MDP}_{i,j} + 0.46}
$$

where $\chi^0$ and $\eta$ are atom parameters, $\mathrm{MDP}$ is the minimal distance path, i.e., the
minimal number of bonds between two atoms, and $K$ is 0.556, 0.778, 1.000, 1.053, 1.087, and 1.091 for
paths of 1 to 6 bonds and 1.091 for longer paths.

The molecular graph must be connected.

(veem)=
#### VEEM

Valence Electrons Equalization Method ([Wu 2011](https://doi.org/10.1088/1674-0068/24/01/31-39))
calculates atomic charges from the numbers of valence electrons. First, the equalized electronegativity
is calculated for the whole molecule:

$$
\chi_{ve} = \frac{\sum_i \chi_i N_{ve,i}}{\sum_i N_{ve,i}}
$$

where $\chi_i$ is the electronegativity of the isolated atom $i$ and $N_{ve,i}$ is its number of valence
electrons. The charge of atom $i$ is:

$$
q_i = N_{ve,i}\, \frac{\chi_{ve} - \chi_i}{\chi_{ve}}
$$

The molecule must be neutral.

### Geometry-dependent methods

These methods also need 3D coordinates and give one set of charges per conformer.

(eem)=
#### EEM

Electronegativity Equalization Method ([Mortier 1986](https://doi.org/10.1021/ja00275a013)) is based on
Sanderson's principle: the electronegativity of each atom gets equalized when atoms bond to form a
molecule. The electronegativity of an atom in a molecule is expressed as:

$$
\chi_i = A_i + B_i q_i + \kappa \sum_{j \neq i} \frac{q_j}{R_{i,j}}
$$

where $A$ and $B$ are atom parameters and $\kappa$ is a common parameter. Setting every $\chi_i$ equal
to $\bar\chi$, together with the total charge:

$$
Q = \sum_i q_i
$$

yields a system of $N + 1$ linear equations:

$$
\begin{bmatrix}
B_1 & \kappa R_{1,2}^{-1} & \cdots & \kappa R_{1,N}^{-1} & 1 \\
\kappa R_{2,1}^{-1} & B_2 & \cdots & \kappa R_{2,N}^{-1} & 1 \\
\vdots & \vdots & \ddots & \vdots & \vdots \\
\kappa R_{N,1}^{-1} & \kappa R_{N,2}^{-1} & \cdots & B_N & 1 \\
1 & 1 & \cdots & 1 & 0
\end{bmatrix}
\begin{bmatrix} q_1 \\ q_2 \\ \vdots \\ q_N \\ -\bar\chi \end{bmatrix}
=
\begin{bmatrix} -A_1 \\ -A_2 \\ \vdots \\ -A_N \\ Q \end{bmatrix}
$$

EEM supports charged molecules and [fixed ions](fixed-ions.md#point-source-terms).

(abeem)=
#### ABEEM

Atom-Bond Electronegativity Equalization Method ([Yang 1997](https://doi.org/10.1021/jp9711048)) extends
[EEM](#eem) to also include bond electronegativities in the equalization. The electronegativity of atom
$i$ is expressed as:

$$
\chi_i = a_i + b_i q_i + c_i \sum_{i-j} q_{i-j} + k \sum_{j \neq i} \frac{q_j}{R_{i,j}}
  + k \sum_{k-l \neq i-j} \frac{q_{k-l}}{R_{i,k-l}}
$$

where $a$, $b$, and $c$ are atom parameters, $k$ is a common parameter, and $q_{i-j}$ is the charge on
the bond $i-j$. The electronegativity of a bond $i-j$ has the form:

$$
\chi_{i-j} = A_{i-j} + B_{i-j} q_{i-j} + C_{i-j} q_i + D_{i-j} q_j
  + k \sum_{l \neq i,j} \frac{q_l}{R_{i-j,l}} + k \sum_{k-l \neq i-j} \frac{q_{k-l}}{R_{i-j,k-l}}
$$

where $A$, $B$, $C$, and $D$ are bond parameters; $D$ applies to the atom with the higher atomic number.
Atomic and bond charges together sum to $Q$.

Each bond charge is placed at the point on the bond whose distances from the two nuclei are proportional
to their covalent radii. Half of each bond charge is then added to each of its two atoms, giving the
final atomic charges.

(sfkeem)=
#### SFKEEM

Selfconsistent Functional Kernel Equalized Electronegativity Method
([Chaves 2006](https://doi.org/10.1021/ci050505e)) builds on [EEM](#eem) with a different hardness matrix:

$$
\chi_i = A_i + 2 B_i q_i + \sum_{j \neq i} 2\sqrt{B_i B_j}\, \operatorname{sech}(\sigma R_{i,j})\, q_j
$$

where $A$ and $B$ are atom parameters, $\sigma$ is a common parameter, and $\operatorname{sech}$ is the
hyperbolic secant.

(qeq)=
#### QEq

Charge Equilibration ([Rappé 1991](https://doi.org/10.1021/j100161a070)) is similar to [EEM](#eem):

$$
\chi_i = \chi_i^0 + J_i q_i + \sum_{j \neq i} J_{i,j} q_j
$$

where $\chi^0$ and $J$ are atom parameters in eV. Originally, the hardness of hydrogen depended on its
charge, which made QEq iterative; ChargeFW uses fixed parameters and solves the system once.

In the original publication, the Coulomb term $J_{i,j}$ was obtained from *ab initio* calculations.
Several empirical terms were developed to replace it with a simple expression. The default
(DasGupta-Huzinaga) follows Oda and Hirono (*J. Mol. Struct. THEOCHEM* 634, 159, 2003):

$$
J_{i,j} = \frac{k}{R_{i,j} + k \big/ \left(\frac{1}{2}(J_i + J_j)\, e^{0.4 R_{i,j} / a_0}\right)}
$$

where $k = 14.4$ eV·Å and $a_0 = 0.529177$ Å. The Louwen-Vogt term uses the geometric mean of Njo, Fan,
and van de Graaf (*J. Mol. Catal. A* 134, 79, 1998):

$$
J_{i,j} = \frac{k}{\left[(k / \sqrt{J_i J_j})^3 + R_{i,j}^3\right]^{1/3}}
$$

The Nishimoto-Mataga, Nishimoto-Mataga-Weiss, Ohno, and Ohno-Klopman terms are also available.

`overlap_term` selects the term (default `DasGupta-Huzinaga`).

(smpqeq)=
#### SMP/QEq

Self-Consistent Charge Equilibration Method ([Zhang 2009](https://doi.org/10.1021/jp8063273)) builds on
[QEq](#qeq) with a hardness that is a function of charge, which makes the scheme iterative:

$$
\chi_i = A_i + 2\lambda(q_i)\, q_i + \sum_{j \neq i} J_{i,j} q_j
$$

where

$$
\lambda(q_i) = B_i + \frac{C_i}{4} q_i + \frac{D_i}{12} q_i^2
$$

and

$$
J_{i,j} = \frac{14.4}{\left[\left(14.4 / 2\sqrt{B_i B_j}\right)^3 + R_{i,j}^3\right]^{1/3}}
$$

where $A$, $B$, $C$, and $D$ are atom parameters. Starting from zero charges, $\lambda$ is evaluated with
the charges of the previous step until no charge changes by more than $10^{-10}$. A calculation that has
not converged after 100 steps fails.

(eqeq)=
#### EQeq

Extended Charge Equilibration Method ([Wilmer 2012](https://doi.org/10.1021/jz3008485)) builds on
[QEq](#qeq). ChargeFW implements the non-periodic case with neutral charge centers:

$$
\chi_i = \chi_i^0 + J_i q_i + \frac{\lambda k}{2} \sum_{j \neq i} q_j \left(\frac{1}{R_{i,j}} + O_{i,j}\right)
$$

where

$$
\chi_i^0 = \frac{\mathrm{IP}_i + \mathrm{EA}_i}{2}, \qquad J_i = \mathrm{IP}_i - \mathrm{EA}_i
$$

IP and EA are the ionization potential and electron affinity in eV, with EA $= -2.0$ eV for hydrogen,
$k = 14.4$ eV·Å, and $\lambda = 1.2$. The overlap term is:

$$
O_{i,j} = e^{-a^2 R_{i,j}^2} \left(2a - a^2 R_{i,j} - \frac{1}{R_{i,j}}\right), \qquad
a = \frac{\sqrt{J_i J_j}}{k}
$$

(eqeqc)=
#### EQeq+C

Bond-order-corrected Extended Charge Equilibration Method
([Martin-Noble 2015](https://doi.org/10.1021/acs.jctc.5b00037)) follows the same procedure as
[EQeq](#eqeq) and then adds corrections to the charges:

$$
q_i = q_i^0 + \sum_{j \neq i} T_{i,j} B_{i,j}
$$

where $q_i^0$ is the EQeq charge and

$$
T_{i,j} = D_i - D_j, \qquad B_{i,j} = e^{-\alpha (R_{i,j} - r_i - r_j)}
$$

where $D$ is an atom parameter, $\alpha$ is a common parameter, and $r$ is the covalent radius. The sum
runs over all atom pairs; the bonds of the input are not used.

(sqe)=
#### SQE

Split-charge Equilibration ([Nistor 2006](https://doi.org/10.1063/1.2346671)) is based on electronegativity
equalization. Unlike EEM or QEq, it works in the bond domain with split charges, i.e., charges located on
the bonds. The split charges $\mathbf{p}$ and the atomic charges $\mathbf{q}$ are:

$$
\left(\mathbf{T} \mathbf{H} \mathbf{T}^{\mathsf T} + \operatorname{diag}(\boldsymbol{\kappa})\right) \mathbf{p}
= -\mathbf{T} \boldsymbol{\chi}, \qquad
\mathbf{q} = \mathbf{T}^{\mathsf T} \mathbf{p}
$$

where $\mathbf{T}$ is the incidence matrix describing the molecular topology,
$\operatorname{diag}(\boldsymbol{\kappa})$ is the diagonal matrix of bond hardnesses,
$\boldsymbol{\chi}$ is the vector of atomic electronegativities, and $\mathbf{H}$ is the hardness matrix
with the atomic hardnesses on the diagonal and the interactions:

$$
H_{i,j} = \frac{\operatorname{erf}\left(R_{i,j} \big/ \sqrt{2 w_i^2 + 2 w_j^2}\right)}{R_{i,j}}
$$

where $w$ is an atom parameter for the charge width.

Every connected component must be formally neutral; use [SQE+q0](#sqeq0) for charged molecules. SQE
supports [fixed ions](fixed-ions.md#point-source-terms).

(sqeq0)=
#### SQE+q0

SQE+q0 ([Verstraelen 2012](https://doi.org/10.1021/ct200512e)) extends [SQE](#sqe) by using the formal
charges $\mathbf{q}_0$ as initial seeds:

$$
\left(\mathbf{T} \mathbf{H} \mathbf{T}^{\mathsf T} + \operatorname{diag}(\boldsymbol{\kappa})\right) \mathbf{p}
= -\mathbf{T} \left(\boldsymbol{\chi} + \mathbf{H} \mathbf{q}_0 - \boldsymbol{\eta} * \mathbf{q}_0\right),
\qquad
\mathbf{q} = \mathbf{T}^{\mathsf T} \mathbf{p} + \mathbf{q}_0
$$

where $\boldsymbol{\eta}$ is the vector of atomic hardnesses and $*$ is the element-wise product. Each
connected component keeps the sum of its formal charges.

SQE+q0 supports charged molecules and [fixed ions](fixed-ions.md#point-source-terms).

(sqeqp)=
#### SQE+qp

SQE+qp ([Schindler 2021](https://doi.org/10.1186/s13321-021-00528-w)) replaces the formal charges with
parameterized ones, $\mathbf{q}_p$. To preserve the total charge of the molecule, they are normalized
first:

$$
\mathbf{q}_{pn} = \mathbf{q}_p - \frac{1}{N}\left(\mathbf{1}^{\mathsf T} \mathbf{q}_p - Q\right)
$$

where $\mathbf{1}^{\mathsf T}$ is a row vector of ones. The rest of the procedure is the same as in
[SQE+q0](#sqeq0), with $\mathbf{q}_{pn}$ in place of $\mathbf{q}_0$. The normalization applies to the
whole molecule, so the parameterized charges of one connected component need not sum to its formal
charge.

SQE+qp supports charged molecules and [fixed ions](fixed-ions.md#point-source-terms).

## Calculation cost

The **Time** column above and the memory reported by `chargefw methods` describe how the cost of one
full calculation grows with the number of atoms $n$ and bonds $m$. Memory counts the method's working
storage, not the molecule itself. These terms show how a method scales; they do not predict seconds or
bytes.

Before calculating, ChargeFW also estimates the working memory and the number of operations of each
plan for your molecules. Dense linear algebra is counted in floating-point operations and other loops
per atom, bond, or atom pair; iterative methods are counted at their maximum number of iterations.
Use the estimates to compare plans, not to predict run time;
[resource estimates](large-structures.md#threads-and-estimates) describes what they include.
