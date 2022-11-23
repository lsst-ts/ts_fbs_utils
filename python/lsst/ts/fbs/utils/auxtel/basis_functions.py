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

__all__ = [
    "get_basis_functions_image_survey",
    "get_basis_functions_cwfs_survey",
    "get_basis_functions_spectroscopic_survey",
]

import typing
from rubin_sim.scheduler import basis_functions


def get_basis_functions_image_survey(
    ra: float,
    nside: int,
    note: str,
    ha_limits: typing.List[typing.Tuple[float, float]],
    wind_speed_maximum: float,
    nobs_reference: int,
    nobs_survey: int,
    note_interest: str,
    filter_names: list,
    gap_min: float,
) -> typing.List[basis_functions.Base_basis_function]:
    """Get the basis functions for the image survey.

    Parameters
    ----------
    ra : `float`
        Right Ascention for the observavtions, in degrees.
    nside : `int`
        The nside value for the healpix grid.
    note : `str`
        A identifier string to be attached to the survey observations.
    ha_limits : 'list` of `float`
        A two-element list with the hour angle limits, in hours.
    wind_speed_maximu : `float`
        Maximum wind speed tolerated for the observations of the survey,
        in m/s.
    nobs_reference : `int`
        Reference number of observations.
    nobs_survey : `int`
        Total number of observations expected for the survey.
    note_interest : `str`
        A substring that maps to surveys to be accounted for against the
        reference number of observations.
    filter_names : list [str]
         List of filter names that need be observed before activating.
    gap_min : `float`
        Gap between subsequent observations, in minutes.

    Returns
    -------
    `list` of `basis_functions.Base_basis_function`
    """

    sun_alt_limit = -12.0

    return [
        basis_functions.Not_twilight_basis_function(sun_alt_limit=sun_alt_limit),
        basis_functions.Hour_Angle_limit_basis_function(RA=ra, ha_limits=ha_limits),
        basis_functions.Slewtime_basis_function(nside=nside, filtername="g"),
        basis_functions.Slewtime_basis_function(nside=nside, filtername="r"),
        basis_functions.Slewtime_basis_function(nside=nside, filtername="i"),
        basis_functions.Moon_avoidance_basis_function(nside=nside),
        basis_functions.Zenith_shadow_mask_basis_function(
            min_alt=20.0, max_alt=85.0, nside=nside
        ),
        basis_functions.VisitGap(note=note, filter_names=filter_names, gap_min=gap_min),
        basis_functions.AvoidDirectWind(
            wind_speed_maximum=wind_speed_maximum, nside=nside
        ),
        basis_functions.BalanceVisits(
            nobs_reference=nobs_reference, note_survey=note, note_interest=note_interest
        ),
        basis_functions.RewardNObsSequence(
            n_obs_survey=nobs_survey,
            note_survey=note,
            nside=nside,
        ),
    ]


def get_basis_functions_cwfs_survey(
    nside: int,
    note: str,
    time_gap_min: float,
    wind_speed_maximum: float,
) -> typing.List[basis_functions.Base_basis_function]:
    """Get the basis functions for the CWFS survey.

    This is a background survey that will activate at specific points in time
    to make sure the telescope optics are aligned.

    Parameters
    ----------
    nside : `int`
        The nside value for the healpix grid.
    note : `str`
        A identifier string to be attached to the survey observations.
    time_gap_min : `float`
        The gap between observations of this survey, in minutes.
    wind_speed_maximu : `float`
        Maximum wind speed tolerated for the observations of the survey,
        in m/s.

    Returns
    -------
    `list` of `basis_functions.Base_basis_function`
    """
    sun_alt_limit = -12.0

    return [
        basis_functions.Not_twilight_basis_function(sun_alt_limit=sun_alt_limit),
        basis_functions.Slewtime_basis_function(nside=nside, filtername="g"),
        basis_functions.Slewtime_basis_function(nside=nside, filtername="r"),
        basis_functions.Slewtime_basis_function(nside=nside, filtername="i"),
        basis_functions.Moon_avoidance_basis_function(nside=nside),
        basis_functions.Zenith_shadow_mask_basis_function(
            min_alt=28.0, max_alt=85.0, nside=nside
        ),
        basis_functions.VisitGap(note=note, gap_min=time_gap_min),
        basis_functions.AvoidDirectWind(
            wind_speed_maximum=wind_speed_maximum, nside=nside
        ),
    ]


def get_basis_functions_spectroscopic_survey(
    ra: float,
    nside: int,
    note: str,
    ha_limits: typing.List[typing.Tuple[float, float]],
    wind_speed_maximum: float,
    gap_min: float,
    moon_distance: float,
    nobs_reference: int,
    note_interest: str,
) -> typing.List[basis_functions.Base_basis_function]:
    """Get basis functions for spectroscopic survey.

    Parameters
    ----------
    ra : `float`
        RA of the target, in degrees.
    nside : `int`
        Healpix maps resolution.
    note : `str`
        Survey note.
    ha_limits : `list` of `tuple` of (`float`, `float`)
        Hour angle limits, in hours.
    wind_speed_maximum : `float`
        Maximum wind speed, in m/s.
    gap_min : `float`
        Gap between subsequent observations, in minutes.
    nobs_reference : `int`
        Reference number of observations.
    note_interest : `str`
        A substring that maps to surveys to be accounted for against the
        reference number of observations.

    Returns
    -------
    list of basis_functions.Base_basis_function
        List of basis functions.
    """
    sun_alt_limit = -12.0

    return [
        basis_functions.Not_twilight_basis_function(sun_alt_limit=sun_alt_limit),
        basis_functions.Hour_Angle_limit_basis_function(RA=ra, ha_limits=ha_limits),
        basis_functions.M5_diff_basis_function(nside=nside),
        basis_functions.Slewtime_basis_function(nside=nside, filtername="g"),
        basis_functions.Slewtime_basis_function(nside=nside, filtername="r"),
        basis_functions.Slewtime_basis_function(nside=nside, filtername="i"),
        basis_functions.Moon_avoidance_basis_function(
            nside=nside, moon_distance=moon_distance
        ),
        basis_functions.Zenith_shadow_mask_basis_function(
            min_alt=28.0, max_alt=85.0, nside=nside
        ),
        basis_functions.VisitGap(note=note, gap_min=gap_min),
        basis_functions.AvoidDirectWind(
            wind_speed_maximum=wind_speed_maximum, nside=nside
        ),
        basis_functions.BalanceVisits(
            nobs_reference=nobs_reference, note_survey=note, note_interest=note_interest
        ),
    ]
