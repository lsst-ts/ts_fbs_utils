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

from rubin_scheduler.scheduler import basis_functions


def get_basis_functions_image_survey(
    ra: float,
    nside: int,
    note: str,
    ha_limits: typing.List[typing.Tuple[float, float]] | None,
    wind_speed_maximum: float,
    nobs_reference: int,
    nobs_survey: int,
    note_interest: str | None,
    filter_names: list,
    gap_min: float,
    additional_notes: list[tuple[str, int]] | None = None,
    avoid_wind: bool = True,
    include_slew: bool = True,
    sun_alt_limit: float = -12,
) -> typing.List[basis_functions.BaseBasisFunction]:
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
    additional_notes : `list` [ `tuple` [`str`, `int`] ], optional
        Optional list of additional notes to configure a second
        RewardNObsSequence basis function to reward completing a set of tiles.
        The first element should be the substring that defines the tile set,
        and the second should be the number of tiles.
    avoid_wind : `bool`, optional
        If True, add the wind avoidance basis function.
        If False, drop basis function entirely.
        Makes use align with the spectroscopic survey.
    include_slew : `bool`, optional
        Include slewtime basis functions (or not).
    sun_alt_limit : `float`, optional
        Sun altitude limit for the survey (degrees).
        Sun must be below this limit for survey to be feasible.

    Returns
    -------
    `list` of `basis_functions.BaseBasisFunction`
    """

    bfs = [
        basis_functions.SunAltLimitBasisFunction(alt_limit=sun_alt_limit),
        basis_functions.MoonAvoidanceBasisFunction(nside=nside),
        basis_functions.AltAzShadowMaskBasisFunction(
            min_alt=26.0, max_alt=85.0, nside=nside, shadow_minutes=0.0, pad=0.0
        ),
        # Note that band_names should include ONLY the bands in use, as
        # the VisitGap will not trigger until all bands are satisfied
        basis_functions.VisitGap(note=note, filter_names=filter_names, gap_min=gap_min),
    ]

    if avoid_wind:
        bfs.append(
            basis_functions.AvoidDirectWind(
                wind_speed_maximum=wind_speed_maximum, nside=nside
            )
        )

    if ha_limits is not None:
        bfs.append(
            basis_functions.HourAngleLimitBasisFunction(RA=ra, ha_limits=ha_limits)
        )

    if include_slew:
        bfs.extend(
            [
                basis_functions.SlewtimeBasisFunction(nside=nside, filtername="g"),
                basis_functions.SlewtimeBasisFunction(nside=nside, filtername="r"),
                basis_functions.SlewtimeBasisFunction(nside=nside, filtername="i"),
            ]
        )

    # The basis functions below are really only relevant for image
    # surveys based on using Tiles for dithering.
    if nobs_survey > 0:
        bfs.append(
            basis_functions.RewardNObsSequence(
                n_obs_survey=nobs_survey,
                note_survey=note,
                nside=nside,
            )
        )

    if note_interest is not None:
        bfs.append(
            basis_functions.BalanceVisits(
                nobs_reference=nobs_reference,
                note_survey=note,
                note_interest=note_interest,
            )
        )

    if additional_notes is not None:
        for additional_note, additional_nobs in additional_notes:
            bfs.append(
                basis_functions.RewardNObsSequence(
                    n_obs_survey=additional_nobs,
                    note_survey=additional_note,
                    nside=nside,
                ),
            )

    return bfs


def get_basis_functions_cwfs_survey(
    nside: int,
    note: str,
    time_gap_min: float,
    wind_speed_maximum: float,
    sun_alt_limit: float = -7,
) -> typing.List[basis_functions.BaseBasisFunction]:
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
    wind_speed_maximum : `float`
        Maximum wind speed tolerated for the observations of the survey,
        in m/s.
    sun_alt_limit : `float`, optional
        Sun altitude limit for the survey (degrees).
        Sun must be below this limit for survey to be feasible.

    Returns
    -------
    `list` of `basis_functions.BaseBasisFunction`
    """

    return [
        basis_functions.M5DiffBasisFunction(nside=nside),
        basis_functions.SlewtimeBasisFunction(nside=nside, filtername="g"),
        basis_functions.SlewtimeBasisFunction(nside=nside, filtername="r"),
        basis_functions.SlewtimeBasisFunction(nside=nside, filtername="i"),
        basis_functions.MoonAvoidanceBasisFunction(nside=nside),
        basis_functions.AltAzShadowMaskBasisFunction(
            min_alt=26.0, max_alt=85.0, nside=nside
        ),
        basis_functions.VisitGap(note=note, gap_min=time_gap_min),
        basis_functions.AvoidDirectWind(
            wind_speed_maximum=wind_speed_maximum, nside=nside
        ),
        basis_functions.SunAltLimitBasisFunction(alt_limit=sun_alt_limit),
    ]


def get_basis_functions_spectroscopic_survey(
    ra: float,
    nside: int,
    note: str,
    ha_limits: typing.List[typing.Tuple[float, float]],
    avoid_wind: bool,
    wind_speed_maximum: float,
    gap_min: float,
    moon_distance: float,
    nobs_reference: int,
    note_interest: str,
    include_slew: bool = True,
    sun_alt_limit: float = -10,
) -> typing.List[basis_functions.BaseBasisFunction]:
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
    avoid_wind : `bool`
        If True, include AvoidDirectWind basis function
    include_slew: `bool`
        If True, include slewtime basis functions
    wind_speed_maximum : `float`
        Maximum wind speed, in m/s.
    gap_min : `float`
        Gap between subsequent observations, in minutes.
    nobs_reference : `int`
        Reference number of observations.
    note_interest : `str`
        A substring that maps to surveys to be accounted for against the
        reference number of observations.
    include_slew : `bool`, optional
        Include slewtime basis functions (or not).
    sun_alt_limit : `float`, optional
        Sun altitude limit for the survey (degrees).
        Sun must be below this limit for survey to be feasible.

    Returns
    -------
    list of basis_functions.BaseBasisFunctio
        List of basis functions.
    """

    bfs = [
        basis_functions.SunAltLimitBasisFunction(alt_limit=sun_alt_limit),
        basis_functions.MoonAvoidanceBasisFunction(
            nside=nside, moon_distance=moon_distance
        ),
        basis_functions.AltAzShadowMaskBasisFunction(
            min_alt=26.0, max_alt=85.0, shadow_minutes=0.0, pad=0.0, nside=nside
        ),
        basis_functions.VisitGap(note=note, gap_min=gap_min),
    ]

    if ha_limits is not None:
        bfs.append(
            basis_functions.HourAngleLimitBasisFunction(RA=ra, ha_limits=ha_limits)
        )

    if avoid_wind:
        bfs.append(
            basis_functions.AvoidDirectWind(
                wind_speed_maximum=wind_speed_maximum, nside=nside
            )
        )

    if include_slew:
        bfs.extend(
            [
                basis_functions.SlewtimeBasisFunction(nside=nside, filtername="g"),
                basis_functions.SlewtimeBasisFunction(nside=nside, filtername="r"),
                basis_functions.SlewtimeBasisFunction(nside=nside, filtername="i"),
            ]
        )
    return bfs
