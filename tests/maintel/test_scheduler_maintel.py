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

from lsst.ts.fbs.utils import Tiles
from lsst.ts.fbs.utils.maintel.make_scheduler import MakeScheduler, SurveyType


class TestMakeScheduler(unittest.TestCase):
    make_scheduler: MakeScheduler

    @classmethod
    def setUpClass(cls) -> None:
        cls.make_scheduler = MakeScheduler()

    def test_get_scheduler_sit(self) -> None:
        nside = 32
        wind_speed_maximum = 8.0
        nside, scheduler = self.make_scheduler.get_scheduler(
            nside=nside,
            wind_speed_maximum=wind_speed_maximum,
            survey_type=SurveyType.SIT,
            image_tiles=[
                Tiles(
                    survey_name="BLOCK-306",
                    hour_angle_limit=[(-6.0, 6.0)],
                    reward_value=0.5,
                    filters=["g", "r", "i"],
                    visit_gap=1440.0,
                    exptime=60.0,
                    nexp=2,
                )
            ],
        )
