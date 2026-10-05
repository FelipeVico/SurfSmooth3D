# Source attribution and external dependencies

`LICENSE` and `COPYRIGHT` preserve the notices supplied by fmm3dbie for its
imported code. The latter refers to original upstream paths; import mappings
in `provenance` identify the corresponding locations in SurfSmooth3D.

- The smoother, MATLAB geometry utilities and basis routines were copied from
  the recorded fmm3dbie working tree. RV node/weight tables retain the original
  attribution to Zydrunas Gimbutas and the upstream exception notices.
- The STEP implementation was copied from the current STEP V3 development
  working tree, including the modified and untracked source files recorded in
  the provenance snapshot. Its source notices were retained; no separate
  package-level license file was present in that source working tree.
- FMM3D remains an external dependency. Its supplied license identifies
  Apache-2.0 and per-file exceptions. Its own LICENSE and AUTHORS are retained
  under `provenance/licenses` for the reference installation and copied
  provenance-only helpers.
- STEP's dependency lock records OCCT 7.9.3 as
  `LGPL-2.1-only WITH OCCT-exception-1.0` and Gmsh 4.12.2 as
  `GPL-2.0-or-later WITH Gmsh-exception`. Their complete upstream notices remain
  in the verified dependency source archives/build area. They are external
  build dependencies rather than committed source copies in this repository.

These records describe the supplied source notices; no authorship or license
terms were changed as part of the extraction.
