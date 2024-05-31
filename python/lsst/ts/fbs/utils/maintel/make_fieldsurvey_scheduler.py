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

from rubin_scheduler.scheduler.schedulers import CoreScheduler
from rubin_scheduler.scheduler.surveys import FieldSurvey

from ..data import field_survey_centers


class MakeFieldSurveyScheduler:
    """Class to construct Feature Based Schedulers consisting of an ensemble of field surveys."""

    def __init__(
        self, 
        nside: int = 32,
    ) -> None:
        
        self.nside = 32
        self.surveys = []


    def _load_candidate_fields(self) -> typing.Dict:
        """Docstring."""
        name_fields, ra_fields, dec_fields = zip(*field_survey_centers.fields)
        fields = {}
        for name, ra, dec in zip(name_fields, ra_fields, dec_fields):
            fields[name] = {"RA": ra, "dec": dec}
        return fields

    
    def add_field_surveys(
        self,
        field_names: typing.List[str],
        **kwargs,
    ) -> None:
        """Add a set of field surveys to the scheduler configuration."""

        # Assert that survey_name is not included among the kwargs
        
        basis_functions = []

        fields = self._load_candidate_fields()
        
        for field_name in field_names:
            RA = fields[field_name]["RA"]
            dec = fields[field_name]["dec"]

            # Consider adding flexibility to add a prefix or suffix
            survey_name = field_name
            
            self.surveys.append(
                FieldSurvey(
                    basis_functions,
                    RA,
                    dec,
                    survey_name=survey_name,
                    **kwargs,
                )
            )


    def add_cwfs_survey(self) -> None:
        """Add curvature wavefront sensing survey."""
        pass

    
    def get_scheduler(self) -> typing.Tuple[int, CoreScheduler]:
        """Construct feature based scheduler for ensemble of field surveys."""

        scheduler = CoreScheduler(self.surveys, nside=self.nside)
        
        return self.nside, scheduler
