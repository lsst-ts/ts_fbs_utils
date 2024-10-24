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

__all__ = [
    "get_basis_functions_star_tracker_survey",
    "get_basis_functions_blob_survey",
    "get_basis_functions_ddf_survey",
    "get_basis_functions_anytime_survey",
]

from rubin_scheduler.scheduler import basis_functions
from rubin_scheduler.scheduler.utils import EuclidOverlapFootprint


def get_basis_functions_star_tracker_survey(
    ra: float,
    nside: int,
    note: str,
    ha_limits: list[tuple[float, float]],
    wind_speed_maximum: float,
    nobs_reference: int,
    nobs_survey: int,
    note_interest: str,
    filter_names: list,
    gap_min: float,
) -> list[basis_functions.BaseBasisFunction]:
    """Get the basis functions for the image survey.

    Parameters
    ----------
    ra : `float`
        Right Ascension for the observations, in degrees.
    nside : `int`
        The nside value for the healpix grid.
    note : `str`
        A identifier string to be attached to the survey observations.
    ha_limits : 'list` of `float`
        A two-element list with the hour angle limits, in hours.
    wind_speed_maximum : `float`
        Maximum wind speed tolerated for the observations of the survey,
        in m/s.
    nobs_reference : `int`
        Reference number of observations.
    nobs_survey : `int`
        Total number of observations expected for the survey.
    note_interest : `str`
        A substring that maps to surveys to be accounted for against the
        reference number of observations.
    filter_names : `list` [ `str` ]
         List of filter names that need be observed before activating.
    gap_min : `float`
        Gap between subsequent observations, in minutes.

    Returns
    -------
    `list` of `basis_functions.BaseBasisFunction`
    """

    sun_alt_limit = -12.0

    return [
        basis_functions.NotTwilightBasisFunction(sun_alt_limit=sun_alt_limit),
        basis_functions.HourAngleLimitBasisFunction(RA=ra, ha_limits=ha_limits),
        basis_functions.SlewtimeBasisFunction(nside=nside, filtername="g"),
        basis_functions.MoonAvoidanceBasisFunction(nside=nside),
        basis_functions.AltAzShadowMaskBasisFunction(
            min_alt=26.0, max_alt=85.0, nside=nside
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


def get_basis_functions_blob_survey(
    nside: int,
    wind_speed_maximum: float,
    footprint: object,
) -> list[basis_functions.BaseBasisFunction]:
    """Get the basis functions for the blob survey.

    Parameters
    ----------
    nside : `int`
        The nside value for the healpix grid.
    wind_speed_maximum : `float`
        Maximum wind speed tolerated for the observations of the survey,
        in m/s.
    footprint: `object`
        Footprint object to generate the blob.

    Returns
    -------
    `list` of `basis_functions.BaseBasisFunction`
    """

    sun_alt_limit = -12.0

    return [
        basis_functions.NotTwilightBasisFunction(sun_alt_limit=sun_alt_limit),
        basis_functions.M5DiffBasisFunction(filtername="r", nside=nside),
        basis_functions.FootprintBasisFunction(
            filtername="r", footprint=footprint, nside=nside
        ),
        basis_functions.MoonAvoidanceBasisFunction(nside=nside),
        basis_functions.AltAzShadowMaskBasisFunction(
            min_alt=26.0, max_alt=85.0, nside=nside
        ),
        basis_functions.AvoidDirectWind(
            wind_speed_maximum=wind_speed_maximum, nside=nside
        ),
        basis_functions.SlewtimeBasisFunction(nside=nside, filtername="r"),
        basis_functions.VisitRepeatBasisFunction(nside=nside),
    ]


def get_basis_functions_ddf_survey(
    nside: int,
    survey_name: str,
    ra: float,
    ha_limits: list[tuple[float, float]],
    wind_speed_maximum: float,
    gap_min: float,
) -> list[basis_functions.BaseBasisFunction]:
    """Get the basis functions for the DDF survey.

    Parameters
    ----------
    nside : `int`
        The nside value for the healpix grid.
    survey_name : `str`
        Name of survey.
    ra : `float`
        Right ascension of target field in degrees.
    ha_limits : 'list` of `float`
        A two-element list with the hour angle limits, in hours.
    wind_speed_maximum : `float`
        Maximum wind speed tolerated for the observations of the survey,
        in m/s.
    gap_min : `float`
        Gap between subsequent observations, in minutes.

    Returns
    -------
    `list` of `basis_functions.BaseBasisFunction`
    """
    sun_alt_limit = -12.0

    return [
        basis_functions.NotTwilightBasisFunction(sun_alt_limit=sun_alt_limit),
        basis_functions.HourAngleLimitBasisFunction(RA=ra, ha_limits=ha_limits),
        basis_functions.AltAzShadowMaskBasisFunction(
            min_alt=26.0, max_alt=85.0, nside=nside
        ),
        basis_functions.AvoidDirectWind(
            wind_speed_maximum=wind_speed_maximum, nside=nside
        ),
        basis_functions.VisitGap(note=survey_name, gap_min=gap_min),
    ]


def get_basis_functions_anytime_survey(
    nside: int,
) -> list[basis_functions.BaseBasisFunction]:
    """Get basis functions for the anytime survey.

    Parameters
    ----------
    nside : `int`
        The healpix map resolution.

    Returns
    -------
    `list`[ `basis_functions.BaseBasisFunction` ]
        List of basis functions.
    """
    sky = EuclidOverlapFootprint()
    footprints, labels = sky.return_maps()
    target_map = footprints["r"]

    bfs = [
        basis_functions.HaMaskBasisFunction(
            ha_min=-1.5,
            ha_max=1.5,
            nside=nside,
        ),
        basis_functions.AltAzShadowMaskBasisFunction(
            min_alt=40.0,
            max_alt=82.0,
            nside=nside,
        ),
        basis_functions.SlewtimeBasisFunction(filtername="r", nside=nside),
        basis_functions.TargetMapBasisFunction(target_map=target_map),
    ]

    return bfs
