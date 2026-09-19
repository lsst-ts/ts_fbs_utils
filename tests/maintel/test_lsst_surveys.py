# This file is part of ts_fbs_utils.
#
# Developed for the Vera Rubin Observatory Telescope and Site System.
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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import unittest

import healpy as hp
import lsst.ts.fbs.utils.maintel.lsst_footprints as lsst_footprints
import lsst.ts.fbs.utils.maintel.lsst_surveys as lsst_surveys
import lsst.ts.fbs.utils.maintel.roman_surveys as roman_surveys
import lsst.ts.fbs.utils.maintel.too_surveys as too_surveys
import numpy as np
import rubin_scheduler.scheduler.basis_functions as basis_functions
from rubin_scheduler.scheduler.surveys import LongGapSurvey
from rubin_scheduler.scheduler.utils import CurrentAreaMap


class Test_lsst_Surveys(unittest.TestCase):
    def setUp(self) -> None:
        # Generate footprints over the sky
        self.nside = 32
        self.footprints, self.template_fp, footprint_mask = (
            lsst_footprints.get_footprints(self.nside)
        )

    def test_footprints(self) -> None:
        footprints, template_fp, footprint_mask = lsst_footprints.get_footprints(
            self.nside
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
