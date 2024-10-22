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

import typing

import astropy.units as u
import numpy as np
from rubin_scheduler.scheduler.detailers import BaseDetailer
from rubin_scheduler.scheduler.surveys import BaseSurvey, BlobSurvey, FieldSurvey
from rubin_scheduler.scheduler.utils import ObservationArray
from rubin_scheduler.utils import ddf_locations

from ..target import Target
from .basis_functions import (
    get_basis_functions_blob_survey,
    get_basis_functions_ddf_survey,
    get_basis_functions_star_tracker_survey,
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

    basis_functions = get_basis_functions_star_tracker_survey(
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
    )

    sequence = [ObservationArray(n=1) for i in range(len(target.filters))]

    for filter_obs, observation in zip(target.filters, sequence):
        observation["RA"] = target.ra.to(u.rad).value
        observation["dec"] = target.dec.to(u.rad).value
        observation["filter"] = filter_obs
        observation["exptime"] = target.exptime
        observation["nexp"] = target.nexp
        observation["scheduler_note"] = f"{target.survey_name}:{target.target_name}"
        observation["target_name"] = target.target_name
        observation["science_program"] = target.survey_name

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
        survey_name=target.survey_name,
        target_name=target.target_name,
        science_program=target.survey_name,
        reward_value=target.reward_value,
        nside=nside,
        nexp=target.nexp,
        detailers=survey_detailers,
    )

    image_survey.basis_weights *= target.reward_value

    return image_survey


def generate_blob_survey(
    nside: int,
    wind_speed_maximum: float,
    footprints: object,
    filter_names: str,
    survey_name: str,
) -> BaseSurvey:
    """Generate blob survey.

    Parameters
    ----------
    nside : `int`
        Healpix map resolution.
    wind_speed_maximum : `float`
        Wind speed limit, in m/s.
    footprints : `object`
        Footprint object to generate the blob.
    filter_names : `list` of `str`
        List of filters to add to the blob survey.
    survey_name : `str`
        Name of survey.

    Returns
    -------
    blob_survey : `BlobSurvey`
        Blob survey.
    """

    basis_functions = get_basis_functions_blob_survey(
        nside=nside, wind_speed_maximum=wind_speed_maximum, footprint=footprints
    )

    basis_weights = np.ones(len(basis_functions))

    blob_survey = BlobSurvey(
        basis_functions,
        basis_weights,
        filtername1=filter_names,
        survey_name="Area",
        science_program=survey_name,
    )

    return blob_survey


def generate_ddf_surveys(
    nside: int,
    wind_speed_maximum: float,
    gap_min: float,
    survey_base_name: str,
) -> BaseSurvey:
    """Generate DDF survey.

    Parameters
    ----------
    nside : `int`
        Healpix map resolution.
    wind_speed_maximum : `float`
        Wind speed limit, in m/s.
    gap_min : `float`
        Gap between subsequent observations, in minutes.
    survey_base_name: `str`
        Name of survey.

    Returns
    -------
    ddf_survey : `DeepDrillingSurvey`
        DDF survey as a DeepDrillingSurvey object.
    """

    locations = ddf_locations()

    # ELAIS S1
    target_field = "DD_ELAISS1"
    ra = locations["ELAISS1"][0]
    dec = locations["ELAISS1"][1]
    ha_limits = [(0.0, 4.5), (19.5, 24.0)]
    basis_functions = get_basis_functions_ddf_survey(
        nside=nside,
        survey_name="DD",
        ra=ra,
        ha_limits=ha_limits,
        wind_speed_maximum=wind_speed_maximum,
        gap_min=gap_min,
    )

    ddf_survey_1 = FieldSurvey(
        basis_functions,
        ra,
        dec,
        sequence="r",
        nvisits=dict(r=20),
        exptimes=dict(r=30),
        survey_name="DD",
        target_name=target_field,
        science_program=survey_base_name,
        nside=nside,
        nexps=dict(r=2),
        detailers=None,
    )

    # Galactic Bulge
    target_field = "DD_GALACTIC_CENTER"
    ra = 270.33
    dec = -27.5
    ha_limits = [(0.0, 4.5), (19.5, 24.0)]
    basis_functions = get_basis_functions_ddf_survey(
        nside=nside,
        survey_name="DD",
        ra=ra,
        ha_limits=ha_limits,
        wind_speed_maximum=wind_speed_maximum,
        gap_min=gap_min,
    )

    ddf_survey_2 = FieldSurvey(
        basis_functions,
        ra,
        dec,
        sequence="r",
        nvisits=dict(r=20),
        exptimes=dict(r=30),
        survey_name="DD",
        target_name=target_field,
        science_program=survey_base_name,
        nside=nside,
        nexps=dict(r=2),
        detailers=None,
    )

    return [ddf_survey_1, ddf_survey_2]
