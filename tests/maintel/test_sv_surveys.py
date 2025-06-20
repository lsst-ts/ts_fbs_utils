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

import os
import unittest

import lsst.ts.fbs.utils.maintel.sv_config as sv_config
import lsst.ts.fbs.utils.maintel.sv_surveys as sv_surveys
import rubin_scheduler.scheduler.basis_functions as basis_functions
from lsst.ts.fbs.utils import get_data_dir
from rubin_scheduler.scheduler.surveys import LongGapSurvey


class Test_SV_Surveys(unittest.TestCase):
    def setUp(self) -> None:
        self.survey_info = sv_config.survey_footprint(survey_start_mjd=60846.5)

    def test_safety_masks(self) -> None:
        nside = 32
        masks = sv_surveys.safety_masks(nside=nside)
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

    def test_generate_blobs(self) -> None:
        surveys = sv_surveys.generate_blobs(footprints=self.survey_info["Footprints"])
        assert len(surveys) > 0

    def test_generate_twi_blobs(self) -> None:
        surveys = sv_surveys.generate_twi_blobs(
            footprints=self.survey_info["Footprints"]
        )
        assert len(surveys) > 0

    def test_gen_template_surveys(self) -> None:
        surveys = sv_surveys.gen_template_surveys(
            footprints=self.survey_info["Footprints"]
        )
        assert len(surveys) > 0

    def test_gen_greedy_surveys(self) -> None:
        surveys = sv_surveys.gen_greedy_surveys(
            footprints=self.survey_info["Footprints"]
        )
        assert len(surveys) > 0

    def test_generate_blobs(self) -> None:
        surveys = sv_surveys.generate_blobs(footprints=self.survey_info["Footprints"])
        assert len(surveys) > 0

    def test_gen_ddf_surveys(self) -> None:
        surveys = sv_surveys.gen_ddf_surveys(
            ddf_config_file=os.path.join(get_data_dir(), "ddf_sv.dat")
        )
        assert len(surveys) > 0

    def test_gen_long_gaps_survey(self) -> None:
        surveys = sv_surveys.gen_long_gaps_survey(
            footprints=self.survey_info["Footprints"]
        )
        assert len(surveys) > 0
        assert isinstance(surveys[0], LongGapSurvey)

    def test_gen_lvk_templates(self) -> None:
        surveys = sv_surveys.gen_lvk_templates(
            footprints_hp=self.survey_info["fp_array"]
        )
        assert len(surveys) > 0
