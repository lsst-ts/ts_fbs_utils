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
import typing
import unittest

from lsst.ts.fbs.utils.auxtel.basis_functions import (
    get_basis_functions_cwfs_survey,
    get_basis_functions_image_survey,
    get_basis_functions_spectroscopic_survey,
)
from rubin_scheduler.scheduler.basis_functions import BaseBasisFunction
from rubin_scheduler.scheduler.model_observatory import ModelObservatory
from rubin_scheduler.scheduler.utils import empty_observation


class TestBasisFunctions(unittest.TestCase):
    def setUp(self) -> None:
        self.obs = empty_observation(1)
        mo = ModelObservatory()
        self.conditions = mo.return_conditions()

    def check_bf_list(self, bf_list: typing.List[BaseBasisFunction]) -> None:
        """Check the basis functions in a list
        have working APIs
        """
        for bf in bf_list:
            bf.add_observation(self.obs)
            feas = bf.check_feasibility(self.conditions)
            reward = bf(self.conditions)
            assert feas is not None
            assert reward is not None

    def test_get_basis_functions_image_survey(self) -> None:
        basis_functions = get_basis_functions_image_survey(
            ra=0.0,
            nside=32,
            note="unit_test_123",
            ha_limits=[(-6.0, 6.0)],
            wind_speed_maximum=6.0,
            nobs_reference=3,
            nobs_survey=30,
            note_interest="unit_test",
            filter_names=["g", "r", "i"],
            gap_min=144.0,
            additional_notes=[["unit_test", 32]],
        )

        assert len(basis_functions) == 12

        self.check_bf_list(basis_functions)

    def test_get_basis_functions_cwfs_survey(self) -> None:
        basis_functions = get_basis_functions_cwfs_survey(
            nside=32,
            note="cwfs",
            time_gap_min=144.0,
            wind_speed_maximum=7.0,
        )

        assert len(basis_functions) == 9
        self.check_bf_list(basis_functions)

    def test_get_basis_functions_spectroscopic_survey(self) -> None:
        basis_functions = get_basis_functions_spectroscopic_survey(
            ra=0.0,
            nside=32,
            note="unit_test_321",
            ha_limits=[(-6.0, 6.0)],
            wind_speed_maximum=4.0,
            avoid_wind=True,
            gap_min=3600.0,
            note_interest="spec",
            moon_distance=30.0,
            nobs_reference=3,
        )

        assert len(basis_functions) == 10
        self.check_bf_list(basis_functions)

    def test_get_basis_functions_spectroscopic_survey_no_wind(self) -> None:
        basis_functions = get_basis_functions_spectroscopic_survey(
            ra=0.0,
            nside=32,
            note="unit_test_321",
            ha_limits=[(-6.0, 6.0)],
            wind_speed_maximum=4.0,
            avoid_wind=False,
            gap_min=3600.0,
            note_interest="spec",
            moon_distance=30.0,
            nobs_reference=3,
        )

        assert len(basis_functions) == 9
        self.check_bf_list(basis_functions)


if __name__ == "__main__":
    unittest.main()
