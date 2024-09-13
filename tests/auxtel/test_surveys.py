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

from astropy import units
from astropy.coordinates import Angle
from lsst.ts.fbs.utils import Target
from lsst.ts.fbs.utils.auxtel.surveys import generate_cwfs_survey, generate_image_survey
from rubin_scheduler.scheduler.model_observatory import ModelObservatory
from rubin_scheduler.scheduler.surveys import BaseSurvey
from rubin_scheduler.scheduler.utils import empty_observation


class TestSurveys(unittest.TestCase):
    def setUp(self) -> None:
        self.obs = empty_observation(1)
        mo = ModelObservatory()
        self.conditions = mo.return_conditions()

    def check_survey(self, survey: BaseSurvey) -> None:
        """Check a survey API works"""

        survey.add_observation(self.obs)
        reward = survey.calc_reward_function(self.conditions)
        obs_out = survey.generate_observations(self.conditions)

        assert reward is not None
        assert obs_out is not None

    def test_generate_image_survey(self) -> None:
        target = Target(
            target_name="unit_test_target",
            survey_name="unit_test_survey",
            ra=Angle("00:00:00", unit=units.hourangle),
            dec=Angle("-30:00:00", unit=units.degree),
            hour_angle_limit=[(-6.0, 6.0)],
            reward_value=1.0,
            filters=["r"],
            visit_gap=144.0,
            exptime=30.0,
            nexp=2,
        )
        survey = generate_image_survey(
            nside=32,
            target=target,
            nfields=64,
            wind_speed_maximum=5.0,
            survey_detailers=[],
        )

        assert survey.survey_name == target.survey_name
        self.check_survey(survey)

    def test_generate_cwfs_survey(self) -> None:
        survey = generate_cwfs_survey(
            nside=32,
            time_gap_min=120.0,
            wind_speed_maximum=8.0,
            cwfs_block_name="cwfs",
        )

        assert survey.survey_name == "cwfs"

        self.check_survey(survey)


if __name__ == "__main__":
    unittest.main()
