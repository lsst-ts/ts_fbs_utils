.. _Version_History:

===============
Version History
===============

v0.4.1
------

* Add target field to auxtel_tiles.txt data file.
* Edit version history to match tag-released version history.  
* Update .github/workflows/lint.yaml python version to 3.10

v0.4.0
------

* Update rubin_sim dependencies to be consistent with v1.0+.
* Update conda build.

v0.3.1
------

* Update zenith_shadow_mask basis function min_alt for all surveys.

v0.3.0
------

* Add optional filter_list keyword to VisitGap imaging survey.
* Update conda build recipe.
* In ```test_basis_functions.py``` update unit tests with filter_list keyword.

v0.2.1
------

* Reduce auxtel imaging survey target fields.

v0.2.0
------

* Add optional moon_distance attribute to target class.
* Add BalanceVisits basis function to spectroscopic survey class.

v0.1.1
------

* Fix text encoding in some targets in auxtel_tiles data file.
* Add support for building conda packages with both python 3.8 and 3.10.

v0.1.0
------

* Initial version of Feature Based Scheduler Utility package focused on AuxTel.
