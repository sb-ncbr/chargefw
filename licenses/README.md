# Third-party notices

ChargeFW is distributed under the MIT License in [`LICENSE`](../LICENSE). Its binaries include code
from the components below; each directory holds their license texts. Which components a binary
contains depends on how it was built: Python wheels, for example, do not contain CLI11, and native
builds without Python bindings do not contain nanobind. The versions are those that ChargeFW pins for
downloaded builds; builds against installed packages use the installed versions.

| Component | Version | License | Source | Files |
| --- | --- | --- | --- | --- |
| Gemmi | 0.7.4 | MPL-2.0 | <https://github.com/project-gemmi/gemmi> | [`gemmi/`](gemmi/) |
| PEGTL (bundled with Gemmi) | as bundled | MIT | <https://github.com/taocpp/PEGTL> | [`gemmi/third_party/`](gemmi/third_party/) |
| fast_float (bundled with Gemmi) | as bundled | Apache-2.0, MIT, or BSL-1.0 | <https://github.com/fastfloat/fast_float> | [`gemmi/third_party/`](gemmi/third_party/) |
| half (bundled with Gemmi) | 2.2.1 | MIT | <https://sourceforge.net/projects/half/> | [`gemmi/third_party/`](gemmi/third_party/) |
| pocketfft (bundled with Gemmi) | as bundled | BSD-3-Clause | <https://github.com/mreineck/pocketfft> | [`gemmi/third_party/`](gemmi/third_party/) |
| tinydir (bundled with Gemmi) | as bundled | BSD-2-Clause | <https://github.com/cxong/tinydir> | [`gemmi/third_party/`](gemmi/third_party/) |
| sajson (bundled with Gemmi) | as bundled | MIT | <https://github.com/chadaustin/sajson> | [`gemmi/third_party/`](gemmi/third_party/) |
| stb_sprintf (bundled with Gemmi) | 1.10 | MIT or public domain | <https://github.com/nothings/stb> | [`gemmi/third_party/`](gemmi/third_party/) |
| QCProt (adapted in Gemmi) | 1.5 | BSD-3-Clause | <https://theobald.brandeis.edu/qcp/> | [`gemmi/third_party/`](gemmi/third_party/) |
| ksw2 (adapted in Gemmi) | as adapted | MIT | <https://github.com/lh3/ksw2> | [`gemmi/third_party/`](gemmi/third_party/) |
| Larch Cromer-Liberman code (adapted in Gemmi) | as adapted | MIT¹ | <https://github.com/xraypy/xraylarch> | [`gemmi/third_party/`](gemmi/third_party/) |
| oneTBB | 2023.1.0 | Apache-2.0 | <https://github.com/uxlfoundation/oneTBB> | [`oneTBB/`](oneTBB/) |
| Eigen | 5.0.1 | MPL-2.0 | <https://gitlab.com/libeigen/eigen> | [`Eigen/`](Eigen/) |
| nanoflann | 2.0.1 | BSD-2-Clause | <https://github.com/jlblancoc/nanoflann> | [`nanoflann/`](nanoflann/) |
| nlohmann/json | 3.12.0 | MIT | <https://github.com/nlohmann/json> | [`nlohmann_json/`](nlohmann_json/) |
| CLI11 | 2.7.2 | BSD-3-Clause | <https://github.com/CLIUtils/CLI11> | [`CLI11/`](CLI11/) |
| nanobind | 2.11 or newer | BSD-3-Clause | <https://github.com/wjakob/nanobind> | [`nanobind/`](nanobind/) |
| robin-map (bundled with nanobind) | as bundled | MIT | <https://github.com/Tessil/robin-map> | [`nanobind/`](nanobind/) |

¹ Gemmi credits the adapted Larch code as 2-clause BSD; Larch is now distributed under the MIT License
included here.

nlohmann/json's [`THIRD-PARTY-NOTICES.txt`](nlohmann_json/THIRD-PARTY-NOTICES.txt) lists the
authors of code it incorporates. Eigen's [`COPYING.README`](Eigen/COPYING.README) describes which of
its license texts apply to which files. oneTBB's [`third-party-programs.txt`](oneTBB/third-party-programs.txt) lists the components it
includes.
