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
import yaml
from rubin_scheduler.scheduler.detailers import BaseDetailer, TrackingInfoDetailer
from rubin_scheduler.scheduler.surveys import BaseSurvey, FieldSurvey, GreedySurvey
from rubin_scheduler.scheduler.utils import ObservationArray

from ..target import Target
from ..utils import get_data_dir
from .basis_functions import (
    get_basis_functions_cwfs_survey,
    get_basis_functions_image_survey,
    get_basis_functions_spectroscopic_survey,
)


def get_auxtel_targets(infile: str | None = None) -> dict:
    """Load potential targets for auxtel observations.

    Targets can be split into different categories for convenience.
    The targets in a particular category should be destined for the same
    'kind' of survey (imaging survey or spectroscopy survey), so that they
    can be configured similarly.

    After target_dict is returned, remove or include desired targets.

    Returns
    -------
    target_pointings : `dict`
        Dictionary of candidate target information.
        The minimum necessary is ra, dec, and block.
        Additional options include exptime, nexp, visit_gap, and priority.
    """
    if infile is None:
        infile = get_data_dir() / "auxtel_targets.yaml"
    with open(infile) as stream:
        target_pointings = yaml.safe_load(stream)
    return target_pointings


def generate_image_survey_from_tiles(
    nside: int,
    target: Target,
    wind_speed_maximum: float,
    nfields: int,
    survey_detailers: typing.List[BaseDetailer],
) -> BaseSurvey:
    """Generate image survey, where original Targets coming from tiles

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
    # One difference in a series of image surveys using tiles
    # vs a single image survey using dithering is in the basis functions
    basis_functions = get_basis_functions_image_survey(
        ra=target.ra.to(unit=u.deg).value,
        nside=nside,
        note=target.survey_name,
        note_interest=target.science_program,
        ha_limits=target.hour_angle_limit,
        wind_speed_maximum=wind_speed_maximum,
        nobs_reference=nfields,
        nobs_survey=len(target.filters),
        filter_names=target.filters,
        gap_min=target.visit_gap,
        additional_notes=[(target.target_name.split("_")[0], 32)],
    )

    # Set up a sequence of nexp exposures per filter
    sequence = [ObservationArray(n=1) for i in range(len(target.filters))]
    for filter_obs, observation in zip(target.filters, sequence):
        observation["RA"] = target.ra.to(u.rad).value
        observation["dec"] = target.dec.to(u.rad).value
        observation["filter"] = filter_obs
        observation["exptime"] = target.exptime
        observation["nexp"] = target.nexp
        observation["scheduler_note"] = target.survey_name

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
        target_name=f"{target.target_name}",
        science_program=f"{target.science_program}",
        scheduler_note=f"{target.survey_name}",
        nside=nside,
        detailers=survey_detailers,
    )
    # Weight all of the basis functions up or down to weight overall survey
    image_survey.basis_weights *= target.reward_value

    return image_survey


def generate_image_survey_from_target(
    nside: int,
    target: Target,
    wind_speed_maximum: float,
    survey_detailers: typing.List[BaseDetailer],
    avoid_wind: bool = True,
    include_slew: bool = True,
    sun_alt_limit: float = -12,
) -> BaseSurvey:
    """Generate image survey, for single Target with dithers

    Parameters
    ----------
    nside : `int`
        Healpix map resolution.
    target : `Target`
        Target for the image survey.
    wind_speed_maximum : `float`
        Wind speed limit, in m/s.
    survey_detailers : `list` of `detailers.BaseDetailer`
        List of survey detailers.
    avoid_wind : `bool`, optional
        If True, add the wind avoidance basis function.
        If False, drop basis function entirely.
        Makes use align with the spectroscopic survey.
    include_slew : `bool`, optional
        If True, include slewtime basis functions.
    sun_alt_limit : `float`, optional
        Sun altitude limit for the survey (degrees).
        Sun must be below this limit for survey to be feasible.

    Returns
    -------
    image_survey : `FieldSurvey`
        Image survey.
    """

    basis_functions = get_basis_functions_image_survey(
        ra=target.ra.to(unit=u.deg).value,
        nside=nside,
        note=target.target_name,
        note_interest=None,
        ha_limits=target.hour_angle_limit,
        wind_speed_maximum=wind_speed_maximum,
        nobs_reference=0,
        nobs_survey=0,
        filter_names=target.filters,
        gap_min=target.visit_gap,
        additional_notes=None,
        include_slew=include_slew,
        sun_alt_limit=sun_alt_limit,
    )

    # Set up a sequence of nexp exposures per filter
    sequence = [ObservationArray(n=target.nexp) for i in range(len(target.filters))]
    for filter_obs, observation in zip(target.filters, sequence):
        observation["RA"] = target.ra.to(u.rad).value
        observation["dec"] = target.dec.to(u.rad).value
        observation["filter"] = filter_obs
        observation["exptime"] = target.exptime
        observation["nexp"] = 1
        observation["scheduler_note"] = target.survey_name

    image_survey = FieldSurvey(
        basis_functions=basis_functions,
        RA=np.array(
            [
                target.ra.to(u.degree).value,
            ]
        ),
        dec=np.array(
            [
                target.dec.to(u.degree).value,
            ]
        ),
        sequence=sequence,
        survey_name=f"{target.survey_name}",
        target_name=f"{target.target_name}",
        science_program=f"{target.science_program}",
        scheduler_note=f"{target.survey_name}",
        nside=nside,
        detailers=survey_detailers,
    )
    # Weight all of the basis functions up or down to weight overall survey
    image_survey.basis_weights *= target.reward_value

    return image_survey


def generate_cwfs_survey(
    nside: int,
    time_gap_min: float,
    wind_speed_maximum: float,
    cwfs_block_name: str,
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
    cwfs_block_name : `str`
        Name of the cwfs block survey.


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
        np.ones_like(basis_functions) * 1.0,
        nside=nside,
        survey_name="cwfs",
        science_program=cwfs_block_name,
        nexp=4,
    )


def generate_spectroscopic_survey(
    nside: int,
    target: Target,
    avoid_wind: bool,
    wind_speed_maximum: float,
    nfields: int,
    survey_detailers: typing.List[BaseDetailer],
    include_slew: bool = True,
    sun_alt_limit: float = -10,
) -> BaseSurvey:
    """Generate Spectroscopic Survey.

    Parameters
    ----------
    nside : `int`
        Healpix map resolution.
    target : `Target`
        Target for the image survey.
    avoid_wind : `bool`
        Include AvoidDirectWind basis function.
    wind_speed_maximum : `float`
        Maximum wind speed (in m/s).
    survey_detailers : `list` of `detailers.BaseDetailer`
        List of survey detailers.
    include_slew : `bool`, optional
        If True, include slewtime basis functions.
    sun_alt_limit : `float`, optional
        Sun altitude limit for the survey (degrees).
        Sun must be below this limit for survey to be feasible.

    Returns
    -------
    spectroscopic_survey : `BaseSurvey`
        Spectroscopic survey.
    """
    basis_functions = get_basis_functions_spectroscopic_survey(
        ra=target.ra.to(u.deg).value,
        nside=nside,
        note=target.target_name,
        note_interest=target.survey_name,
        ha_limits=target.hour_angle_limit,
        avoid_wind=avoid_wind,
        wind_speed_maximum=wind_speed_maximum,
        gap_min=target.visit_gap,
        moon_distance=target.moon_distance,
        nobs_reference=nfields,
        include_slew=include_slew,
        sun_alt_limit=sun_alt_limit,
    )

    observation = ObservationArray(n=target.nexp)
    observation["RA"] = target.ra.to(u.rad).value
    observation["dec"] = target.dec.to(u.rad).value
    observation["filter"] = "r"
    observation["exptime"] = target.exptime
    observation["nexp"] = 1
    observation["scheduler_note"] = target.survey_name
    observation = [observation]

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
        sequence=observation,
        survey_name=f"{target.survey_name}",
        target_name=f"{target.target_name}",
        scheduler_note=f"{target.survey_name}",
        science_program=f"{target.science_program}",
        nside=nside,
        detailers=survey_detailers,
    )
    # Weight all of the basis functions up or down to weight overall survey
    spectroscopic_survey.basis_weights *= target.reward_value

    return spectroscopic_survey
