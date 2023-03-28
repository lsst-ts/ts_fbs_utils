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

__all__ = ["MakeScheduler"]

import enum
import typing
import dataclasses

from astropy import units
from astropy.coordinates import Angle

from rubin_sim.scheduler.schedulers import CoreScheduler
from rubin_sim.scheduler.detailers import BaseDetailer
from rubin_sim.scheduler.surveys import BaseSurvey

from .. import get_auxtel_tiles, Target, Tiles, AssertSurvey

from .surveys import (
    generate_cwfs_survey,
    generate_image_survey,
    generate_spectroscopic_survey,
)


class SurveyType(enum.IntEnum):
    Spec = enum.auto()
    Image = enum.auto()
    SpecImage = enum.auto()
    ImageSpec = enum.auto()


class MakeScheduler:
    """A class to help construct Feature Based Schedulers."""

    def __init__(self) -> None:
        self.survey_type_assertions: typing.Dict[SurveyType, AssertSurvey] = {
            SurveyType.Spec: self.assert_spec_survey,
            SurveyType.Image: self.assert_image_survey,
            SurveyType.SpecImage: self.assert_spec_image,
            SurveyType.ImageSpec: self.assert_image_spec,
        }

        self.tiles = get_auxtel_tiles()

    def assert_spec_survey(
        self, spec_targets: typing.List[Target], image_tiles: typing.List[Tiles]
    ) -> None:
        """Assert spectroscopic survey.

        Parameters
        ----------
        spec_targets : `list` of `Target`
            List of targets for spectroscopic survey.
        image_targets : `list` of `Target`
            List of targets for background image survey.

        Raises
        ------
        AssertionError:
            If inputs are invalid.
        """
        assert len(spec_targets) > 0 and len(image_tiles) == 0

    def assert_image_survey(
        self, spec_targets: typing.List[Target], image_tiles: typing.List[Tiles]
    ) -> None:
        """Assert image survey.

        Parameters
        ----------
        spec_targets : `list` of `Target`
            List of targets for spectroscopic survey.
        image_targets : `list` of `Target`
            List of targets for background image survey.

        Raises
        ------
        AssertionError:
            If inputs are invalid.
        """
        assert len(spec_targets) == 0 and len(image_tiles) > 0

    def assert_spec_image(
        self, spec_targets: typing.List[Target], image_tiles: typing.List[Tiles]
    ) -> None:
        """Assert spectroscopic and image survey (spectroscopic with higher
        priority).

        Parameters
        ----------
        spec_targets : `list` of `Target`
            List of targets for spectroscopic survey.
        image_targets : `list` of `Target`
            List of targets for background image survey.

        Raises
        ------
        AssertionError:
            If inputs are invalid.
        """
        assert len(spec_targets) > 0 and len(image_tiles) > 0

    def assert_image_spec(
        self, spec_targets: typing.List[Target], image_tiles: typing.List[Tiles]
    ) -> None:
        """Assert image and spectroscopic survey (image with higher priority).

        Parameters
        ----------
        spec_targets : `list` of `Target`
            List of targets for spectroscopic survey.
        image_targets : `list` of `Target`
            List of targets for background image survey.

        Raises
        ------
        AssertionError:
            If inputs are invalid.
        """
        assert len(spec_targets) > 0 and len(image_tiles) > 0

    def assert_scheduler_inputs(
        self,
        survey_type: SurveyType,
        spec_targets: typing.List[Target],
        image_tiles: typing.List[Tiles],
    ) -> None:
        """Assert that the scheduler inputs are correct.

        Parameters
        ----------
        survey_type : `SurveyType`
            Enumeration with the type of survey. This defines the order of the
            surveys in the scheduler.
        spec_targets : `list` of `Target`
            List of targets for spectroscopic survey.
        image_targets : `list` of `Target`
            List of targets for background image survey.

        Raises
        ------
        AssertionError
            If inputs are invalid.
        """
        self.survey_type_assertions[survey_type](
            spec_targets=spec_targets, image_tiles=image_tiles
        )

    def get_scheduler(
        self,
        nside: int,
        wind_speed_maximum: float,
        survey_type: SurveyType,
        spec_targets: typing.List[Target],
        image_tiles: typing.List[Tiles],
        spec_detailers: typing.List[BaseDetailer],
        image_detailers: typing.List[BaseDetailer],
    ) -> typing.Tuple[int, CoreScheduler]:
        """Construct feature based scheduler for spectroscopic survey with
        image survey in the background (with lower priority).

        Parameters
        ----------
        nside : `int`
            Healpix map resolution.
        wind_speed_maximum : `int`
            Maximum wind speed (in m/s).
        survey_type : `SurveyType`
            Enumeration with the type of survey. This defines the order of the
            surveys in the scheduler.
        spec_targets : `list` of `Target`
            List of targets for spectroscopic survey.
        tiles : `list` of `Tiles`
            List of targets for background image survey.
        spec_detailers : `list` of `BaseDetailer`
            List of Detailers used for spectroscopic survey.
        image_detailers : `list` of `BaseDetailer`
            List of Detailers used for image survey.

        Returns
        -------
        nside : `int`
            Healpix map resolution.
        scheduler : `CoreScheduler`
            Feature based scheduler.
        """

        self.assert_scheduler_inputs(
            survey_type=survey_type,
            spec_targets=spec_targets,
            image_tiles=image_tiles,
        )

        cwfs_survey: typing.List[BaseSurvey] = [
            generate_cwfs_survey(
                nside=nside,
                time_gap_min=120.0,
                wind_speed_maximum=wind_speed_maximum,
            ),
        ]

        spectroscopic_survey: typing.List[BaseSurvey] = []

        for target in spec_targets:
            spectroscopic_survey.append(
                generate_spectroscopic_survey(
                    nside=nside,
                    target=target,
                    wind_speed_maximum=wind_speed_maximum,
                    nfields=len(spec_targets),
                    survey_detailers=spec_detailers,
                )
            )

        # Background imaging survey
        image_target_surveys = self.generate_targets_from_tiles(image_tiles)

        image_survey: typing.List[BaseSurvey] = []

        # Image surveys
        for image_targets in image_target_surveys:
            for target in image_targets:
                image_survey.append(
                    generate_image_survey(
                        nside=nside,
                        target=target,
                        wind_speed_maximum=wind_speed_maximum,
                        nfields=len(image_targets),
                        survey_detailers=image_detailers,
                    )
                )

        first_survey, second_survey = (
            (spectroscopic_survey, image_survey)
            if survey_type in {SurveyType.SpecImage, SurveyType.Spec}
            else (image_survey, spectroscopic_survey)
        )

        surveys = [
            survey
            for survey in [
                cwfs_survey,
                first_survey,
                second_survey,
            ]
            if len(survey) > 0
        ]

        scheduler = CoreScheduler(surveys, nside=nside)

        return nside, scheduler

    def generate_targets_from_tiles(
        self, tiles: typing.List[Tiles]
    ) -> typing.List[typing.List[Target]]:
        """Generate a list of targets from a list of tiles.

        Parameters
        ----------
        image_tiles : `list` of `Tiles`
            List of tiles.

        Returns
        -------
        `list` of `list` of `Target`
            List of targets from the given tiles.
        """
        return [
            [
                Target(
                    target_name=target["Name"],
                    ra=Angle(target["RA"], unit=units.hourangle),
                    dec=Angle(target["Dec"], unit=units.degree),
                    **dataclasses.asdict(tile),
                )
                for target in self.tiles
                if tile.survey_name == target["Survey"]
            ]
            for tile in tiles
        ]
