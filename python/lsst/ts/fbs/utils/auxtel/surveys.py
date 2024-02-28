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

import astropy.units as u
import numpy as np
from rubin_scheduler.scheduler.detailers import BaseDetailer
from rubin_scheduler.scheduler.surveys import BaseSurvey, FieldSurvey, GreedySurvey
from rubin_scheduler.scheduler.utils import empty_observation

from ..target import Target
from .basis_functions import (
    get_basis_functions_cwfs_survey,
    get_basis_functions_image_survey,
    get_basis_functions_spectroscopic_survey,
)


def generate_image_survey(
    nside: int,
    target: Target,
    wind_speed_maximum: float,
    nfields: int,
    survey_detailers: typing.List[BaseDetailer],
) -> BaseSurvey:
    """Generate image survey.

    Parameters
    ----------
    nside : `int`
        Healpix map resolution.
    target : `Target`
        Target for the image survey.
    wind_speed_maximum : `float`
        Wind speed limit, in m/s.
    nfields : `int`
        Number of fields in the general survey.
    survey_detailers : `list` of `detailers.BaseDetailer`
        List of survey detailers.

    Returns
    -------
    image_survey : `FieldSurvey`
        Image survey.
    """

    basis_functions = get_basis_functions_image_survey(
        ra=target.ra.to(unit=u.deg).value,
        nside=nside,
        note=target.target_name,
        note_interest=target.survey_name,
        ha_limits=target.hour_angle_limit,
        wind_speed_maximum=wind_speed_maximum,
        nobs_reference=nfields,
        nobs_survey=len(target.filters),
        filter_names=target.filters,
        gap_min=target.visit_gap,
        additional_notes=[(target.target_name.split("_")[0], 32)],
    )

    sequence = [empty_observation() for i in range(len(target.filters))]

    for filter_obs, observation in zip(target.filters, sequence):
        observation["RA"] = target.ra.to(u.rad).value
        observation["dec"] = target.dec.to(u.rad).value
        observation["filter"] = filter_obs
        observation["exptime"] = target.exptime
        observation["nexp"] = target.nexp
        observation["note"] = f"{target.survey_name}:{target.target_name}"

    image_survey = FieldSurvey(
        basis_functions,
        np.array(
            [
                target.ra.to(u.degree).value,
            ]
        ),
        np.array(
            [
                target.dec.to(u.degree).value,
            ]
        ),
        sequence=sequence,
        survey_name=f"{target.survey_name}",
        reward_value=target.reward_value,
        nside=nside,
        nexp=target.nexp,
        detailers=survey_detailers,
    )

    image_survey.basis_weights *= target.reward_value

    return image_survey


def generate_cwfs_survey(
    nside: int,
    time_gap_min: float,
    wind_speed_maximum: float,
) -> BaseSurvey:
    """Generate Curvature Wavefront Sensing Survey.

    Parameters
    ----------
    nside : `int`
        Healpix map resolution.
    time_gap_min : `float`
        Minimum time between cwfs sequences (in minutes).
    wind_speed_maximum : `float`
        Maximum wind speed (in m/s).

    Returns
    -------
    `GreedySurvey`
        CWFS survey.
    """
    basis_functions = get_basis_functions_cwfs_survey(
        nside=nside,
        note="cwfs",
        time_gap_min=time_gap_min,
        wind_speed_maximum=wind_speed_maximum,
    )

    return GreedySurvey(
        basis_functions,
        np.ones_like(basis_functions) * 10000.0,
        nside=nside,
        survey_name="cwfs",
        nexp=4,
    )


def generate_spectroscopic_survey(
    nside: int,
    target: Target,
    wind_speed_maximum: float,
    nfields: int,
    survey_detailers: typing.List[BaseDetailer],
) -> BaseSurvey:
    """Generate Spectroscopic Survey.

    Parameters
    ----------
    nside : `int`
        Healpix map resolution.
    target : `Target`
        Target for the image survey.
    wind_speed_maximu : `float`
        Maximum wind speed (in m/s).
    survey_detailers : `list` of `detailers.BaseDetailer`
        List of survey detailers.

    Returns
    -------
    spectroscopic_survey : `BaseSurvey`
        Spectroscopic survey.
    """
    basis_functions = get_basis_functions_spectroscopic_survey(
        ra=target.ra.to(u.deg).value,
        nside=nside,
        note=target.target_name,
        note_interest="spec",
        ha_limits=target.hour_angle_limit,
        wind_speed_maximum=wind_speed_maximum,
        gap_min=target.visit_gap,
        moon_distance=target.moon_distance,
        nobs_reference=nfields,
    )

    observation = empty_observation()
    observation["RA"] = target.ra.to(u.rad).value
    observation["dec"] = target.dec.to(u.rad).value
    observation["filter"] = "r"
    observation["exptime"] = target.exptime
    observation["nexp"] = target.nexp
    observation["note"] = f"{target.survey_name}:{target.target_name}"
    sequence = [observation]

    spectroscopic_survey = FieldSurvey(
        basis_functions,
        np.array(
            [
                target.ra.to(u.degree).value,
            ]
        ),
        np.array(
            [
                target.dec.to(u.degree).value,
            ]
        ),
        sequence=sequence,
        survey_name=f"{target.target_name}",  # This is target name on purpose.
        reward_value=target.reward_value,
        nside=nside,
        nexp=target.nexp,
        detailers=survey_detailers,
    )

    spectroscopic_survey.basis_weights *= target.reward_value

    return spectroscopic_survey
