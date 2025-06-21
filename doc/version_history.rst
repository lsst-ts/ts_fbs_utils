.. _Version_History:

===============
Version History
===============

v0.17.1
-------
* Fix some issues with the SV configuration: observation_reason needs underscores and no spaces, and n_obs template should not be 0 for surveys that use rubin_scheduler <= 3.10.0 if they also use the NObsPerYear basis function.

v0.17.0
-------
* Add SV survey configuration support - generation of the footprint, as well as the various tiers of surveys (DDF, long-gaps (triplets), template gathering, standard pairs, twilight pairs, greedy (single) visits, and a final layer of early template gathering). This includes a DDF prescheduling generation file that more likely will long-term live in rubin_scheduler but is currently divergent from the version we have there (on purpose, to suppose "ocean" ddfs).

v0.16.1
-------

* Fix bug in MakeFieldSurveyScheduler to make targets optional.
* Update get_data_dir to allow users setting a altenative data path from an environment variable.
* Added pointing tiles for block 387.

v0.16.0
-------

* Update tools for building fieldsurvey scheduler configurations for LSSTCam
* Move the definition of fieldsurvey pointing centers to ts_config_ocs for more rapid deployment of updates

v0.15.1
-------
* Added SunAltLimitBasisFunctions to cwfs, imaging and spectroscopic surveys, to ensure that the FBS stops requesting targets beyond twilight time. The default values are -10 deg, -12 deg, and -10 deg respectively. Placing these limits lightens the load on the OSs around enabling and disabling the Scheduler near twilight.

v0.15.0
-------

* In auxtel/surveys, add an imaging survey based on a single target, with a single (larger dither) detailer.
* In auxtel/make_scheduler and utils.py - adds a type of imaging survey based on target, in addition to an imaging survey based on tiles (image_target, in addition to image_tiles). 
* In auxtel/basis_functions, remove M5Diff basis functions from imaging and spectroscopy surveys, in order to choose targets based on time when desired rather than when the best m5 is achieved.
* Adds data/auxtel_targets.yaml, to provide a list of targets for auxtel with easy configuration for individual targets.
* In auxtel/surveys, adds a function to read the auxtel_targets.yaml file.
* Extends the target.py data class to include `science_program` (aka json block) explicitly, separate from `survey_name`.

v0.14.3
-------

* In auxtel/basis_functions, override default values for shadow_minutes and pad in AltAzShadowMaskBasisFunction.
* In auxtel/surveys, pass target name as scheduler_note to field survey.

v0.14.2
-------

* Update AuxTel Tiles.
* Update github workflows to pin identify v2.6+

v0.14.1
-------

* Add candidate targets for ComCam science scheduler configurations.

v0.14.0
-------

* Update ``MakeFieldSurveyScheduler`` with ``add_field_altaz_surveys`` method, that will add ``FieldAltAzSurvey`` to the list of surveys.

v0.13.1
-------

* Add candidate targets for ComCam science scheduler configurations.

v0.13.0
-------

* Add utilities for generating ComCam science scheduler configurations.

v0.12.0
-------

* General updates to support migration to rubin-scheduler >2.

v0.11.0
-------

* Replace deprecated ``ZenithShadowMaskBasisFunction`` with ``AltAzShadowMaskBasisFunction``.

v0.10.0
-------

* Update ``auxtel/make_scheduler`` and ``auxtel/surveys`` to allow setting the cwfs survey name.

v0.9.3
------

* In ``data/auxtel_tiles.txt``, remove unused tiles and to rename AUXTEL_PHOTO_IMAGING survey to BLOCK-306.

v0.9.2
------
* In ``data/auxtel_tiles.txt``, expand photo imaging survey.
* Add bool option to toggle AvoidDirectWind basis function for AuxTel spectroscopic surveys.
* Add option to configure cwfs time gap. 

v0.9.1
------
* In ``data/auxtel_tiles.txt``, remove unused target regions.
* Update git lint workflows to call tssw common workflow. 

v0.9.0
------
* Move imports from rubin_sim to rubin_scheduler where applicable.
* Update and extend conda dependencies.

v0.8.0
------

* In ``auxtel/basis_functions.py``, add m5diff basis function to cwfs survey.
* Add MaskAzimuthBasisFucntion to maintel blob and DD surveys.
* Update maintel survey to use ``FieldSurvey`` class instead of ``DeepDrillingSurvey``.
* Add an "anytime" survey to the maintel scheduler to allow using it anytime in the day.

v0.7.3
------

* In ``data/auxtel_tiles.txt``, fixup new target regions for photo-calib survey.

v0.7.2
------

* In ``conda/meta.yaml``, update rubin_sim min_pin build and test requirements.

v0.7.1
------

* In ``data/auxtel_tiles.txt``, add target regions for photo-calib survey.
* Remove unused pytest options from ``pyproject.toml``.
* Update github workflow to skip pre-commit install.

v0.7.0
------

* Add method to configure maintel SIT survey using ``BlobSurvey`` and ``DeepDrillingSurvey``.
* Remove ``conda/conda_build_config.yaml`` file. 

v0.6.1
------

* In ``data/auxtel_tiles.txt``, update target regions for photo-calib survey.

v0.6.0
------

* In ``auxtel/basis_functions``, add second RewardNObs sequence to reward completing a region of tiles
* Update targets in to auxtel_tiles.txt data file for photo-calib survey. 


v0.5.1
------

* In ``auxtel/basis_functions``, remove balanceVisists from spec survey and update unit test. 
* Add pre_commit_config file. 
* Add new targets to auxtel_tiles.txt data file for photo-calib survey. 

v0.5.0
------

* Add option to pass list of detailers to spec and image surveys.
* Update git workflow to run ``pre-commit`` check using ``ts-pre-commit-config``.
* Run ``black`` and ``isort`` in the entire package.
* Update ``.gitignore`` with new ``ts-pre-commit-config``.
* Add sub-module to configure maintel for star tracker survey using ``FieldSurvey``.
* In ``auxtel/basis_functions``, fix typos in docstrings.
* In ``utils``, add new method to read maintel tiles.
* Add tiles for Main Telescope surveys.

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
