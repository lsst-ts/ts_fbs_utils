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

import lsst.ts.fbs.utils.maintel.lsst_surveys as lsst_surveys
import lsst.ts.fbs.utils.maintel.roman_surveys as roman_surveys
import lsst.ts.fbs.utils.maintel.too_surveys as too_surveys
import numpy as np
import rubin_scheduler.scheduler.basis_functions as basis_functions
from rubin_scheduler.scheduler.surveys import LongGapSurvey
from rubin_scheduler.scheduler.utils import CurrentAreaMap, make_rolling_footprints
from rubin_scheduler.site_models import Almanac
from rubin_scheduler.utils import SURVEY_START_MJD


class Test_lsst_Surveys(unittest.TestCase):
    def setUp(self) -> None:
        # Generate footprint over the sky
        self.nside = 32
        sky = CurrentAreaMap(nside=self.nside)
        footprints_hp_array, labels = sky.return_maps()
        # Identify pixels for rolling
        roll_indx = np.where((labels == "lowdust") | (labels == "virgo"))[0]
        roll_footprint = footprints_hp_array["r"] * 0
        roll_footprint[roll_indx] = 1

        footprints_hp = {}
        for key in footprints_hp_array.dtype.names:
            footprints_hp[key] = footprints_hp_array[key]

        # Set up a mask to contain some surveys within this region
        footprint_mask = footprints_hp["r"] * 0
        footprint_mask[np.where(footprints_hp["r"] > 0)] = 1

        # Use Almanac to find the position of the sun at the start of survey
        almanac = Almanac(mjd_start=SURVEY_START_MJD)
        sun_moon_info = almanac.get_sun_moon_positions(SURVEY_START_MJD)
        sun_ra_start = sun_moon_info["sun_RA"].copy()

        # Define the rolling footprint
        footprints = make_rolling_footprints(
            fp_hp=footprints_hp,
            mjd_start=SURVEY_START_MJD,
            sun_ra_start=sun_ra_start,
            nslice=2,
            scale=0.9,
            nside=self.nside,
            wfd_indx=roll_indx,
            order_roll=1,
            n_cycles=3,
            uniform=True,
        )
        self.footprints = footprints

    def test_safety_masks(self) -> None:
        masks = lsst_surveys.safety_masks(nside=self.nside)
        assert len(masks) > 0
        required_masks = [
            basis_functions.MoonAvoidanceBasisFunction,
            basis_functions.AltAzShadowMaskBasisFunction,
            basis_functions.AvoidDirectWind,
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

    def test_generate_short_blobs(self) -> None:
        surveys = lsst_surveys.generate_short_blobs(
            footprints=self.footprints,
            nside=self.nside,
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
