# This file is part of ts_fbs_utils.
#
# Developed for the Vera C. Rubin Observatory Telescope and Site Systems.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

import unittest

import healpy as hp
import lsst.ts.fbs.utils.maintel.lsst_footprints as lsst_footprints
import lsst.ts.fbs.utils.maintel.lsst_surveys as lsst_surveys
import lsst.ts.fbs.utils.maintel.roman_surveys as roman_surveys
import lsst.ts.fbs.utils.maintel.too_surveys as too_surveys
import numpy as np
import rubin_scheduler.scheduler.basis_functions as basis_functions
import rubin_scheduler.scheduler.detailers as detailers
from rubin_scheduler.scheduler.surveys import BaseSurvey, LongGapSurvey
from rubin_scheduler.scheduler.utils import CurrentAreaMap


class Test_lsst_Surveys(unittest.TestCase):
    def setUp(self) -> None:
        # Generate footprints over the sky
        self.nside = 32
        self.footprints, self.template_fp, self.footprint_mask = (
            lsst_footprints.get_footprints(self.nside)
        )

    def test_footprints(self) -> None:
        footprints, template_fp, footprint_mask = (
            self.footprints,
            self.template_fp,
            self.footprint_mask,
        )
        assert len(footprints.footprints[0]) == hp.nside2npix(self.nside)
        # Check that all bands are present in footprints and template_fp
        for band in "ugrizy":
            assert band in self.footprints.bands
            assert band in self.template_fp.bands
        # Check that only the bands in the bandpasses are actually non-zero
        footprints, template_fp, footprint_mask = lsst_footprints.get_footprints(
            self.nside, bandpasses=("g", "r")
        )
        for band in "gr":
            idx = footprints.bands[band]
            assert np.nanmax(footprints.footprints[idx]) > 0
            assert np.nanmax(template_fp.footprints[idx]) > 0
        for band in "uizy":
            idx = footprints.bands[band]
            assert np.nanmax(footprints.footprints[idx]) == 0
            assert np.nanmax(template_fp.footprints[idx]) == 0

    def test_standard_masks(self) -> None:
        masks = lsst_surveys.standard_masks(nside=self.nside)
        assert len(masks) > 0
        required_masks = [
            basis_functions.MoonAvoidanceBasisFunction,
            basis_functions.AltAzShadowMaskBasisFunction,
            basis_functions.MaskDirectWindBasisFunction,
        ]
        # There can be other masks too
        for req_mask in required_masks:
            mask_present = False
            for m in masks:
                if isinstance(m, req_mask):
                    mask_present = True
                    break
            assert mask_present

    def test_standard_bf(self) -> None:
        bfs = lsst_surveys.standard_bf(nside=self.nside, footprints=self.footprints)
        assert len(bfs) > 0

    def test_generate_blobs(self) -> None:
        surveys = lsst_surveys.generate_blobs(
            footprints=self.footprints, nside=self.nside
        )
        assert len(surveys) > 0

    def test_generate_neo_blobs(self) -> None:
        surveys = lsst_surveys.generate_twilight_near_sun(nside=self.nside)
        assert len(surveys) > 0

    def test_gen_template_surveys(self) -> None:
        surveys = lsst_surveys.gen_template_surveys(
            footprints=self.footprints,
            nside=self.nside,
        )
        assert len(surveys) > 0

    def test_gen_greedy_surveys(self) -> None:
        surveys = lsst_surveys.gen_greedy_surveys(footprints=self.footprints)
        assert len(surveys) > 0

    def test_gen_long_gaps_survey(self) -> None:
        surveys = lsst_surveys.gen_long_gaps_survey(
            footprints=self.footprints,
            nside=self.nside,
        )
        assert len(surveys) > 0
        assert isinstance(surveys[0], LongGapSurvey)

    def test_templates(self) -> None:
        surveys = lsst_surveys.gen_template_surveys(
            footprints=self.footprints, nside=self.nside
        )
        assert len(surveys) > 0

    def test_standard_masks_time_limited_shadow(self) -> None:
        masks = lsst_surveys.standard_masks(nside=self.nside)
        assert not any(
            isinstance(m, basis_functions.AltAzShadowTimeLimitedBasisFunction)
            for m in masks
        )

        masks = lsst_surveys.standard_masks(
            nside=self.nside, apply_time_limited_shadow=True
        )
        assert any(
            isinstance(m, basis_functions.AltAzShadowTimeLimitedBasisFunction)
            for m in masks
        )

    def test_non_default_nside(self) -> None:
        nside = 16
        footprints, _, _ = lsst_footprints.get_footprints(nside)

        surveys = (
            lsst_surveys.gen_template_surveys(
                footprints=footprints, nside=nside, band1s=["g"], band2s=["r"]
            )
            + lsst_surveys.generate_blobs(
                footprints=footprints, nside=nside, band1s=["g"], band2s=["r"]
            )
            + lsst_surveys.blob_for_long(footprints=footprints, nside=nside)
        )
        for survey in surveys:
            assert survey.nside == nside, survey.survey_name

    def test_unpaired_bands(self) -> None:
        templates = lsst_surveys.gen_template_surveys(
            footprints=self.footprints,
            nside=self.nside,
            band1s=["u", "g"],
            band2s=None,
        )
        assert [s.survey_name for s in templates] == ["template_u", "template_g"]

        blobs = lsst_surveys.generate_blobs(
            footprints=self.footprints,
            nside=self.nside,
            band1s=["r", "i"],
            band2s=None,
        )
        assert [s.survey_name for s in blobs] == ["pair_33 r", "pair_33 i"]

    def test_mismatched_band_lengths(self) -> None:
        for generator in (
            lsst_surveys.gen_template_surveys,
            lsst_surveys.generate_blobs,
            lsst_surveys.blob_for_long,
        ):
            with self.subTest(generator=generator.__name__):
                with self.assertRaises(ValueError):
                    generator(
                        footprints=self.footprints,
                        nside=self.nside,
                        band1s=["g", "r"],
                        band2s=["i"],
                    )

    def test_extinction_limit(self) -> None:
        def has_cloud_handling(survey: BaseSurvey) -> tuple[bool, bool]:
            has_mask = any(
                isinstance(bf, basis_functions.MaskCloudMapBasisFunction)
                for bf in survey.basis_functions
            )
            has_detailer = any(
                isinstance(det, detailers.ExtinctionLimitDetailer)
                for det in survey.detailers
            )
            return has_mask, has_detailer

        generators = {
            "templates": lambda **kw: lsst_surveys.gen_template_surveys(
                footprints=self.footprints,
                nside=self.nside,
                band1s=["g"],
                band2s=["r"],
                **kw,
            ),
            "blobs": lambda **kw: lsst_surveys.generate_blobs(
                footprints=self.footprints,
                nside=self.nside,
                band1s=["g"],
                band2s=["r"],
                **kw,
            ),
            "greedy": lambda **kw: lsst_surveys.gen_greedy_surveys(
                footprints=self.footprints, nside=self.nside, **kw
            ),
            "twilight_near_sun": lambda **kw: (
                lsst_surveys.generate_twilight_near_sun(nside=self.nside, **kw)
            ),
        }
        for name, generator in generators.items():
            with self.subTest(generator=name):
                for survey in generator(extinction_limit=1.5):
                    assert has_cloud_handling(survey) == (True, True)
                for survey in generator(extinction_limit=None):
                    assert has_cloud_handling(survey) == (False, False)

    def test_template_dither(self) -> None:
        surveys = lsst_surveys.gen_template_surveys(
            footprints=self.footprints, nside=self.nside, band1s=["g"], band2s=["r"]
        )
        assert all(s.dither == "call" for s in surveys)

        blob_survey_params = dict(lsst_surveys.BLOB_SURVEY_PARAMS_DEFAULTS)
        blob_survey_params["dither"] = "night"
        surveys = lsst_surveys.gen_template_surveys(
            footprints=self.footprints,
            nside=self.nside,
            band1s=["g"],
            band2s=["r"],
            blob_survey_params=blob_survey_params,
        )
        assert all(s.dither == "night" for s in surveys)


class Test_too_Surveys(unittest.TestCase):
    def test_too(self) -> None:
        # Generate footprint over the sky
        self.nside = 32
        sky = CurrentAreaMap(nside=self.nside)
        footprints_hp_array, labels = sky.return_maps()

        footprints_hp = {}
        for key in footprints_hp_array.dtype.names:
            footprints_hp[key] = footprints_hp_array[key]

        too_footprint = np.where(footprints_hp["r"] > 0, 1.0, np.nan)

        surveys = too_surveys.gen_too_surveys(
            nside=self.nside, too_footprint=too_footprint
        )
        assert len(surveys) > 0


class Test_roman_surveys(unittest.TestCase):
    def test_roman(self) -> None:
        on_season = roman_surveys.gen_roman_on_season()
        off_season = roman_surveys.gen_roman_off_season()
        assert on_season is not None
        assert off_season is not None
