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

__all__ = ["MakeFieldSurveyScheduler", "get_sv_targets"]

import typing

import yaml
from lsst.ts.fbs.utils import get_pointing_model_grid_data
from rubin_scheduler.scheduler.basis_functions import BaseBasisFunction, VisitGap
from rubin_scheduler.scheduler.detailers import BaseDetailer
from rubin_scheduler.scheduler.schedulers import CoreScheduler
from rubin_scheduler.scheduler.surveys import BaseSurvey, FieldAltAzSurvey, FieldSurvey


class MakeFieldSurveyScheduler:
    """Construct Feature Based Scheduler configuration consisting of an
    ensemble of field surveys, and optionally, other types of surveys.
    """

    def __init__(
        self,
        targets: dict | None = None,
        nside: int = 32,
        ntiers: int = 1,
        band_to_filter: dict | None = None,
    ) -> None:
        self.targets = targets
        self.nside = nside
        self.surveys: typing.List[typing.List[BaseSurvey]] = [
            [] for _ in range(0, ntiers)
        ]
        self.band_to_filter = band_to_filter

    def add_field_surveys(
        self,
        tier: int,
        observation_reason: str,
        science_program: str,
        target_names: typing.List[str],
        survey_name: str | None = None,
        basis_functions: typing.List[BaseBasisFunction] = [],
        sequence: str = "ugrizy",
        nvisits: dict | None = None,
        exptimes: dict | None = None,
        nexps: dict | None = None,
        detailers: typing.List[BaseDetailer] = [],
        **kwargs: typing.Any,
    ) -> None:
        """Add a list of field surveys to the scheduler configuration.

        Parameters
        ----------
        tier: `int`
            Tier index used to control prioritization of surveys.
        json_block: `str`
            JSON BLOCK used to perform the observations.
        observation_reason: `str`
            Purpose of the observation, e.g., "science"
        science_program: `str`
            Name of the science program.
        target_names: `list` of `str`
            List of names to specify the pointing center of each field survey.
        basis_functions: `list` of `BaseBasisFunction`
            Basis functions provided to each field survey.
        """

        for target_name in target_names:
            if self.targets is not None:
                RA = self.targets[target_name]["ra"]
                dec = self.targets[target_name]["dec"]
            else:
                raise ValueError("No list of targets passed.")

            self.surveys[tier].append(
                FieldSurvey(
                    basis_functions,
                    RA,
                    dec,
                    sequence=sequence,
                    nvisits=nvisits,
                    exptimes=exptimes,
                    nexps=nexps,
                    ignore_obs=None,
                    survey_name=survey_name,
                    target_name=target_name,
                    science_program=science_program,
                    observation_reason=observation_reason,
                    scheduler_note=None,
                    nside=self.nside,
                    flush_pad=30.0,
                    detailers=detailers,
                )
            )

    def add_field_altaz_surveys(
        self,
        tier: int,
        observation_reason: str,
        science_program: str,
        survey_name: str | None = None,
        basis_functions: typing.List[BaseBasisFunction] = [],
        sequence: str = "ugrizy",
        nvisits: dict | None = None,
        exptimes: dict | None = None,
        nexps: dict | None = None,
        detailers: typing.List[BaseDetailer] = [],
        **kwargs: typing.Any,
    ) -> None:
        """Add a list of field surveys to the scheduler configuration.

        Parameters
        ----------
        tier: `int`
            Tier index used to control prioritization of surveys.
        json_block: `str`
            JSON BLOCK used to perform the observations.
        observation_reason: `str`
            Purpose of the observation, e.g., "science"
        science_program: `str`
            Name of the science program.
        basis_functions: `list` of `BaseBasisFunction`
            Basis functions provided to each field survey.
        """

        targets = get_pointing_model_grid_data(science_program=science_program)

        for target_name, alt, az in targets:
            self.surveys[tier].append(
                FieldAltAzSurvey(
                    basis_functions=[VisitGap(note=target_name, gap_min=420)]
                    + basis_functions,
                    alt=alt,
                    az=az,
                    sequence=sequence,
                    nvisits=nvisits,
                    exptimes=exptimes,
                    nexps=nexps,
                    ignore_obs=None,
                    survey_name=survey_name,
                    target_name=target_name,
                    science_program=science_program,
                    observation_reason=observation_reason,
                    scheduler_note=target_name,
                    nside=self.nside,
                    flush_pad=30.0,
                    detailers=detailers,
                )
            )

    def add_cwfs_survey(self) -> None:
        """Add curvature wavefront sensing survey."""
        pass

    def get_scheduler(self) -> typing.Tuple[int, CoreScheduler]:
        """Construct feature based scheduler for ensemble of field surveys.

        Returns
        -------
        nside : `int`
            Healpix map resolution.
        scheduler : `CoreScheduler`
            Feature based scheduler.
        """

        scheduler = CoreScheduler(
            self.surveys, nside=self.nside, band_to_filter=self.band_to_filter
        )

        return self.nside, scheduler


def get_sv_targets(target_file: str, exclude: typing.List[str] = []) -> dict:
    """Load candidate target coordinates for field survey science observations.

    Parameters
    ----------
    target_file: `str`
        Filename for YAML file with target names and coordinates.
    exclude : `list[str]`
        List of target names to exclude when loading.

    Returns
    -------
    target_dict : `dict`
        Dictionary of candidate target names and coordinates.
    """
    with open(target_file) as stream:
        targets_dict = yaml.safe_load(stream)

    targets_dict = targets_dict
    return {_: targets_dict[_] for _ in targets_dict if _ not in exclude}
