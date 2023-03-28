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

import dataclasses

from astropy import units
from astropy.coordinates import Angle
from rubin_sim.scheduler.detailers import BaseDetailer
from rubin_sim.scheduler.schedulers import CoreScheduler
from rubin_sim.scheduler.surveys import BaseSurvey

from .. import Target, Tiles, get_maintel_tiles
from .surveys import generate_image_survey


class MakeScheduler:
    """A class to help construct Feature Based Schedulers."""

    def __init__(self) -> None:
        self.tiles = get_maintel_tiles()

    def get_scheduler(
        self,
        nside: int,
        wind_speed_maximum: float,
        image_tiles: list[Tiles],
        survey_detailers: list[BaseDetailer] | None = None,
    ) -> tuple[int, CoreScheduler]:
        """Construct feature based scheduler.

        Parameters
        ----------
        nside : `int`
            Healpix map resolution.
        wind_speed_maximum : `int`
            Maximum wind speed (in m/s).
        image_tiles : `list`[`Tiles`]
            List of tiles.
        survey_details : `list`[`BaseDetailer`], optional
            List of survey detailers (default=None).

        Returns
        -------
        nside : `int`
            Healpix map resolution.
        scheduler : `CoreScheduler`
            Feature based scheduler.
        """

        image_target_surveys = self.generate_targets_from_tiles(image_tiles)

        image_survey: list[BaseSurvey] = []

        # Image surveys
        for image_targets in image_target_surveys:
            for target in image_targets:
                target.visit_gap = 480.0
                image_survey.append(
                    generate_image_survey(
                        nside=nside,
                        target=target,
                        wind_speed_maximum=wind_speed_maximum,
                        nfields=len(image_targets),
                        survey_detailers=survey_detailers
                        if survey_detailers is not None
                        else [],
                    )
                )

        scheduler = CoreScheduler([image_survey], nside=nside)

        return nside, scheduler

    def generate_targets_from_tiles(self, tiles: list[Tiles]) -> list[list[Target]]:
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
