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

import numpy as np
import rubin_scheduler.scheduler.basis_functions as basis_functions
from lsst.ts.fbs.utils.maintel.lsst_footprints import SURVEY_START_MJD
from lsst.ts.fbs.utils.maintel.stability_surveys import (
    gen_az_el_rot_stability_survey,
)
from rubin_scheduler.scheduler.features import Conditions


class TestStabilitySurveys(unittest.TestCase):
    def setUp(self) -> None:
        self.nside = 32
        self.az_values = [0.0, 90.0, 180.0]
        self.el_values = [60.0, 70.0, 80.0]
        self.rot_values = [0.0, 10.0, -10.0]

    def test_gen_az_el_rot_stability_survey(self) -> None:
        surveys = gen_az_el_rot_stability_survey(
            self.az_values,
            self.el_values,
            self.rot_values,
            science_program="BLOCK-T000",
            nside=self.nside,
        )

        assert len(surveys) == len(self.az_values)
        for survey, rot in zip(surveys, self.rot_values):
            np.testing.assert_allclose(
                survey.observations["rotTelPos"], np.radians(rot)
            )

    def test_mismatched_lengths(self) -> None:
        with self.assertRaises(ValueError):
            gen_az_el_rot_stability_survey(
                self.az_values,
                self.el_values[:-1],
                self.rot_values,
                science_program="BLOCK-T000",
                nside=self.nside,
            )

    def test_wind_mask(self) -> None:
        surveys = gen_az_el_rot_stability_survey(
            self.az_values,
            self.el_values,
            self.rot_values,
            science_program="BLOCK-T000",
            nside=self.nside,
        )

        def make_conditions(wind_speed: float) -> Conditions:
            # Conditions caches wind_pressure, so use a fresh object
            # for each wind speed.
            conditions = Conditions(nside=self.nside, mjd=SURVEY_START_MJD + 0.2)
            conditions.wind_speed = wind_speed
            conditions.wind_direction = 0.0
            return conditions

        for survey in surveys:
            wind_masks = [
                bf
                for bf in survey.basis_functions
                if isinstance(bf, basis_functions.MaskDirectWindBasisFunction)
            ]
            assert len(wind_masks) == 1
            wind_mask = wind_masks[0]
            assert wind_mask.wind_speed_maximum is not None

            # Low wind: nothing is masked.
            assert not np.any(np.isnan(wind_mask(make_conditions(5.0))))

            # High wind: the region facing the wind is masked.
            assert np.any(np.isnan(wind_mask(make_conditions(60.0))))


if __name__ == "__main__":
    unittest.main()
