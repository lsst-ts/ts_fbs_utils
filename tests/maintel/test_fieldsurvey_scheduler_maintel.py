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

from lsst.ts.fbs.utils.maintel.make_fieldsurvey_scheduler import (
    MakeFieldSurveyScheduler,
    get_comcam_sv_targets,
)


class TestMakeFieldSurveyScheduler(unittest.TestCase):
    make_scheduler: MakeFieldSurveyScheduler

    @classmethod
    def setUpClass(cls) -> None:
        nside = 32
        ntiers = 1
        cls.make_scheduler = MakeFieldSurveyScheduler(nside, ntiers)

    def test_get_scheduler(self) -> None:
        tier = 0
        observation_reason = "science"
        science_program = "BLOCK-320"
        target_names = get_comcam_sv_targets().keys()
        self.make_scheduler.add_field_surveys(
            tier, observation_reason, science_program, target_names
        )

        nside, scheduler = self.make_scheduler.get_scheduler()


class TestGetComCamSVTargets(unittest.TestCase):

    def test_get_comcam_sv_targets(self) -> None:
        assert len(get_comcam_sv_targets()) > 0

        key = list(get_comcam_sv_targets().keys())[0]
        assert key not in get_comcam_sv_targets(exclude=[key]).keys()
