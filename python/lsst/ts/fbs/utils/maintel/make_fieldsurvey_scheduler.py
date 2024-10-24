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

__all__ = ["MakeFieldSurveyScheduler"]

import typing

from rubin_scheduler.scheduler.basis_functions import BaseBasisFunction
from rubin_scheduler.scheduler.schedulers import CoreScheduler
from rubin_scheduler.scheduler.surveys import BaseSurvey, FieldSurvey

from ..data import field_survey_centers


class MakeFieldSurveyScheduler:
    """Construct Feature Based Scheduler configuration consisting of an
    ensemble of field surveys, and optionally, other types of surveys.
    """

    def __init__(
        self,
        nside: int = 32,
        ntiers: int = 1,
    ) -> None:

        self.nside = 32
        self.surveys: typing.List[typing.List[BaseSurvey]] = [
            [] for _ in range(0, ntiers)
        ]

    def _load_candidate_fields(self) -> typing.Dict:
        """Load pointing center data for field surveys."""
        name_fields, ra_fields, dec_fields = zip(*field_survey_centers.fields)
        fields = {}
        for name, ra, dec in zip(name_fields, ra_fields, dec_fields):
            fields[name] = {"RA": ra, "dec": dec}
        return fields

    def add_field_surveys(
        self,
        tier: int,
        json_block: str,
        observation_reason: str,
        science_program: str,
        field_names: typing.List[str],
        basis_functions: typing.List[BaseBasisFunction] = [],
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
        field_names: `list` of `str`
            List of names to specify the pointing center of each field survey.
        basis_functions: `list` of `BaseBasisFunction`
            Basis functions provided to each field survey.
        """

        if "survey_name" in kwargs:
            raise ValueError("Use program_name rather than survey_name.")

        fields = self._load_candidate_fields()

        for field_name in field_names:
            field_ra_deg = fields[field_name]["RA"]
            field_dec_deg = fields[field_name]["dec"]

            self.surveys[tier].append(
                FieldSurvey(
                    basis_functions,
                    field_ra_deg,
                    field_dec_deg,
                    # sequence=sequence,
                    # nvisits=nvisits,
                    # exptimes=exptimes,
                    # nexps=nexps,
                    ignore_obs=None,
                    survey_name=None,
                    target_name=field_name,
                    science_program=science_program,
                    observation_reason=observation_reason,
                    scheduler_note=None,
                    readtime=2.4,
                    filter_change_time=120.0,
                    nside=self.nside,
                    flush_pad=30.0,
                    # detailers=detailers,
                    # Start old
                    # basis_functions,
                    # field_ra_deg,
                    # field_dec_deg,
                    # json_block=json_block,
                    # observation_reason=observation_reason,
                    # science_program=science_program,
                    # ignore_obs=None,
                    # accept_obs=[field_name],
                    # survey_name=field_name,
                    # scheduler_note=None,
                    # nside=self.nside,
                    # **kwargs,
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

        scheduler = CoreScheduler(self.surveys, nside=self.nside)

        return self.nside, scheduler
