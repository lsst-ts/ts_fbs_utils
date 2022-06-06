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
from astropy import units

from astropy.coordinates import Angle

from lsst.ts.fbs.utils import Target, Tiles
from lsst.ts.fbs.utils.auxtel.make_scheduler import MakeScheduler, SurveyType


class TestMakeScheduler(unittest.TestCase):
    make_scheduler: MakeScheduler
    spec_targets: typing.List[Target]
    image_tiles: typing.List[Tiles]

    @classmethod
    def setUpClass(cls) -> None:
        cls.make_scheduler = MakeScheduler()
        cls.spec_targets = cls.get_spec_targets()
        cls.image_tiles = cls.get_image_tiles()
        return super().setUpClass()

    def test_get_scheduler_spec_image(
        self,
    ) -> None:

        nside, scheduler = self.make_scheduler.get_scheduler(
            nside=32,
            wind_speed_maximum=18.0,
            survey_type=SurveyType.SpecImage,
            spec_targets=self.spec_targets,
            image_tiles=self.image_tiles,
        )

        assert nside == 32
        assert scheduler.nside == nside
        assert len(scheduler.survey_lists) == 3
        assert (
            scheduler.survey_lists[1][0].survey_name == self.spec_targets[0].target_name
        )
        assert (
            scheduler.survey_lists[2][0].survey_name == self.image_tiles[0].survey_name
        )

    def test_get_scheduler_spec_image_fail_no_spec_target(self) -> None:

        with self.assertRaises(AssertionError):
            self.make_scheduler.get_scheduler(
                nside=32,
                wind_speed_maximum=18.0,
                survey_type=SurveyType.SpecImage,
                spec_targets=[],
                image_tiles=self.image_tiles,
            )

    def test_get_scheduler_spec_image_fail_no_image_tiles(self) -> None:

        with self.assertRaises(AssertionError):
            self.make_scheduler.get_scheduler(
                nside=32,
                wind_speed_maximum=18.0,
                survey_type=SurveyType.SpecImage,
                spec_targets=self.spec_targets,
                image_tiles=[],
            )

    def test_get_scheduler_image_spec(self) -> None:

        nside, scheduler = self.make_scheduler.get_scheduler(
            nside=32,
            wind_speed_maximum=18.0,
            survey_type=SurveyType.ImageSpec,
            spec_targets=self.spec_targets,
            image_tiles=self.image_tiles,
        )

        assert nside == 32
        assert scheduler.nside == nside
        assert len(scheduler.survey_lists) == 3
        assert (
            scheduler.survey_lists[1][0].survey_name == self.image_tiles[0].survey_name
        )
        assert (
            scheduler.survey_lists[2][0].survey_name == self.spec_targets[0].target_name
        )

    def test_get_scheduler_image_spec_fail_no_spec_target(self) -> None:

        with self.assertRaises(AssertionError):
            self.make_scheduler.get_scheduler(
                nside=32,
                wind_speed_maximum=18.0,
                survey_type=SurveyType.ImageSpec,
                spec_targets=[],
                image_tiles=self.image_tiles,
            )

    def test_get_scheduler_image_spec_fail_no_image_tiles(self) -> None:

        with self.assertRaises(AssertionError):
            self.make_scheduler.get_scheduler(
                nside=32,
                wind_speed_maximum=18.0,
                survey_type=SurveyType.ImageSpec,
                spec_targets=self.spec_targets,
                image_tiles=[],
            )

    def test_get_scheduler_spec(self) -> None:

        nside, scheduler = self.make_scheduler.get_scheduler(
            nside=32,
            wind_speed_maximum=18.0,
            survey_type=SurveyType.Spec,
            spec_targets=self.spec_targets,
            image_tiles=[],
        )

        assert nside == 32
        assert scheduler.nside == nside
        assert len(scheduler.survey_lists) == 2
        assert (
            scheduler.survey_lists[1][0].survey_name == self.spec_targets[0].target_name
        )

    def test_get_scheduler_spec_fail_with_image_target(self) -> None:

        with self.assertRaises(AssertionError):
            self.make_scheduler.get_scheduler(
                nside=32,
                wind_speed_maximum=18.0,
                survey_type=SurveyType.Spec,
                spec_targets=self.spec_targets,
                image_tiles=self.image_tiles,
            )

    def test_get_scheduler_image(self) -> None:

        nside, scheduler = self.make_scheduler.get_scheduler(
            nside=32,
            wind_speed_maximum=18.0,
            survey_type=SurveyType.Image,
            spec_targets=[],
            image_tiles=self.image_tiles,
        )

        assert nside == 32
        assert scheduler.nside == nside
        assert len(scheduler.survey_lists) == 2
        assert (
            scheduler.survey_lists[1][0].survey_name == self.image_tiles[0].survey_name
        )

    def test_get_scheduler_image_fail_with_spec_target(self) -> None:

        with self.assertRaises(AssertionError):
            self.make_scheduler.get_scheduler(
                nside=32,
                wind_speed_maximum=18.0,
                survey_type=SurveyType.Image,
                spec_targets=self.spec_targets,
                image_tiles=self.image_tiles,
            )

    @staticmethod
    def get_spec_targets() -> typing.List[Target]:
        return [
            Target(
                target_name="TestTarget",
                survey_name="UnitTest",
                ra=Angle("00:00:00.0", unit=units.hourangle),
                dec=Angle("-30:00:00", unit=units.deg),
                hour_angle_limit=[(-6.0, 6.0)],
                reward_value=1.0,
                filters=["r"],
                visit_gap=123.0,
                exptime=120.0,
                nexp=1,
            ),
        ]

    @staticmethod
    def get_image_tiles() -> typing.List[Tiles]:
        return [
            Tiles(
                survey_name="LATISS_POLE",
                hour_angle_limit=[(-6.0, 6.0)],
                reward_value=0.5,
                filters=["g", "r", "i"],
                visit_gap=1440.0,
                exptime=60.0,
                nexp=2,
            )
        ]
