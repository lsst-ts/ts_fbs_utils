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

import numpy as np
import rubin_scheduler.scheduler.basis_functions as bf
import rubin_scheduler.scheduler.detailers as detailers
from rubin_scheduler.utils import DEFAULT_NSIDE

__all__ = [
    "safety_masks",
    "standard_bf"
]

# Set up values to use as kwarg defaults.
NEXP = 1
U_NEXP = 1
EXPTIME = 30.0
U_EXPTIME = 38.0
CAMERA_ROT_LIMITS = (-80.0, 80.0)
SCIENCE_PROGRAM = "BLOCK-365"


def safety_masks(
    nside: int = DEFAULT_NSIDE,
    moon_distance: float = 30,
    wind_speed_maximum: float = 20.0,
    min_alt: float = 20,
    max_alt: float = 86.5,
    min_az: float = 0,
    max_az: float = 360,
    shadow_minutes: float = 70,
    min_az_sunrise: float = 120,
    max_az_sunrise: float = 290,
    time_to_sunrise: float = 3.0,
) -> list[bf.BaseBasisFunction]:
    """Basic safety mask basis functions.

    Avoids the moon, bright planets, high wind, and
    areas on the sky out of bounds, using
    the MoonAvoidanceBasisFunction, PlanetMaskBasisFunction,
    AvoidDirectWindBasisFunction, and the AltAzShadowMaskBasisFunction.
    Adds the default AltAzShadowMaskTimeLimited basis function to avoid
    pointing toward sunrise late in the night during commissioning.

    Parameters
    ----------
    nside : `int` or None
        The healpix nside to use.
        Default of None uses rubin_scheduler.utils.get_default_nside.
    moon_distance : `float`, optional
        Moon avoidance distance, in degrees.
    wind_speed_maximum : `float`, optional
        Wind speed maximum to apply to the wind avoidance basis function,
        in m/s.
    min_alt : `float`, optional
        Minimum altitude (in degrees) to observe.
    max_alt : `float`, optional
        Maximum altitude (in degrees) to observe.
    min_az : `float`, optional
        Minimum azimuth angle (in degrees) to observe.
    max_az : `float`, optional
        Maximum azimuth angle (in degrees) to observe.
    shadow_minutes : `float`, optional
        Avoid inaccessible alt/az regions, as well as parts of the sky
        which will move into those regions within `shadow_minutes` (minutes).
    min_az_sunrise : `float`, optional
        Minimum azimuth angle (in degrees) to observe during time period
        at the end of the night (during time_to_sunrise).
    max_az_sunrise: `float`, optional
        Maximum azimuth angle (in degrees) to observe during the time period
        at the end of the night (during time_to_sunrise).
    time_to_sunrise : `float`, optional
        Hours before daybreak (sun @ alt=0) to start the azimuth avoidance
        mask.

    Returns
    -------
    mask_basis_functions : `list` [`BaseBasisFunction`]
        Mask basis functions should always be used with a weight of 0.
        The masked (np.nan or -np.inf) regions will remain masked,
        but the basis function values won't influence the reward.
    """
    mask_bfs = []
    # Avoid the moon - too close to the moon will trip the REBs
    # Re-evaluate this based upon the October testing
    mask_bfs.append(
        bf.MoonAvoidanceBasisFunction(nside=nside, moon_distance=moon_distance)
    )
    # Avoid bright planets
    mask_bfs.append(bf.PlanetMaskBasisFunction(nside=nside))
    # Avoid the wind
    mask_bfs.append(
        bf.AvoidDirectWind(nside=nside, wind_speed_maximum=wind_speed_maximum)
    )
    # Avoid the alt/az limits - this will pick up limits from the
    # yaml file configurations for the summit as well
    # Re-evaluate the minimum altitude restriction
    mask_bfs.append(
        bf.AltAzShadowMaskBasisFunction(
            nside=nside,
            min_alt=min_alt,
            max_alt=max_alt,
            min_az=min_az,
            max_az=max_az,
            shadow_minutes=shadow_minutes,
        )
    )
    # Only look toward the southeast in the morning,
    # permitting emergency dome closure
    mask_bfs.append(
        bf.AltAzShadowTimeLimitedBasisFunction(
            nside=nside,
            min_alt=min_alt,
            max_alt=max_alt,
            min_az=min_az_sunrise,
            max_az=max_az_sunrise,
            shadow_minutes=shadow_minutes,
            # Time until/after sun_keys in hours
            time_to_sun=time_to_sunrise + shadow_minutes / 60.0,
            # 'sunrise' is 0 degree sunrise
            sun_keys=["sunrise"],
        )
    )
    # We should move this into the basis function itself.
    if shadow_minutes > 40:
        mask_bfs.append(
            bf.AltAzShadowTimeLimitedBasisFunction(
                nside=nside,
                min_alt=min_alt,
                max_alt=max_alt,
                min_az=min_az_sunrise,
                max_az=max_az_sunrise,
                shadow_minutes=shadow_minutes / 2.0,
                # Time until/after sun_keys in hours
                time_to_sun=time_to_sunrise + shadow_minutes / 60.0,
                # 'sunrise' is 0 degree sunrise
                sun_keys=["sunrise"],
            )
        )
    return mask_bfs


def standard_bf(
    nside: int = DEFAULT_NSIDE,
    bandname: str = "g",
    bandname2: str = "i",
    m5_weight: float = 6.0,
    footprint_weight: float = 1.5,
    slewtime_weight: float = 3.0,
    stayband_weight: float = 3.0,
    template_weight: float = 12.0,
    u_template_weight: float = 50.0,
    g_template_weight: float = 50.0,
    footprints: Footprints | None = None,
    fiducial_fwhm: float = 1.3,
    n_obs_template: dict | None = None,
    season: float = 365.25,
    season_start_hour: float = -4.0,
    season_end_hour: float = 2.0,
    strict: bool = True,
) -> list[tuple[bf.BaseBasisFunction, float]]:
    """Generate the standard basis functions that are shared by blob surveys.

    Parameters
    ----------
    nside : `int`
        The HEALpix nside to use. Defaults to DEFAULT_NSIDE
    bandname : `str`
        The band name for the first observation. Default "g".
    bandname2 : `str`
        The band name for the second in the pair (None if unpaired).
        Default "i".
    m5_weight : `float`
        The weight for the 5-sigma depth difference basis function.
    footprint_weight : `float`
        The weight on the survey footprint basis function.
    slewtime_weight : `float`
        The weight on the slewtime basis function.
    stayband_weight : `float`
        The weight on basis function that tries to stay avoid band changes.
    template_weight : `float`
        The weight to place on getting image templates every season.
    u_template_weight : `float`
        The weight to place on getting image templates in u-band. Since there
        are so few u-visits, it can be helpful to turn this up a little
        higher than the standard template_weight kwarg.
    g_template_weight : `float`
        The weight to place on getting image templates in g-band. Since there
        are so few g-visits, it can be helpful to turn this up a
        little higher than the standard template_weight kwarg.
    footprints : `rubin_scheduler.scheduler.utils.Footprints` object
        The desired footprints object. Default of None will work, but is likely
        not desirable.
    fiducial_fwhm : `float`
        The fiducial FWHM for the M5Diff Basis function.
    n_obs_template : `dict`
        The number of observations to take every season in each band.
    season : `float`
        The length of season (i.e., how long before templates expire) (days).
        Default 365.25.
    season_start_hour : `float`
        Hour angle limits to use when gathering templates.
        Default -4 (hours)
    season_end_hour : `float`
       Hour angle limits to use when gathering templates.
       Default +2 (hours)
    strict : `bool`
        If False, use BandChangeBasisFunction which rewards visits in the
        same bandpass as currently in-use.
        If True, use a StrictBandBasisFunction which rewards visits in the
        same bandpass as currently in-use, but also rewards visits in
        different filters if the moon rose/set or twilight ended/started,
        or if there was a large gap in observing.

    Returns
    -------
    basis_functions_weights : `list`
        list of tuple pairs (basis function, weight) that is
        (rubin_scheduler.scheduler.BasisFunction object, float)

    """
    template_weights = {
        "u": u_template_weight,
        "g": g_template_weight,
        "r": template_weight,
        "i": template_weight,
        "z": template_weight,
        "y": template_weight,
    }

    bfs = []

    if bandname2 is not None:
        bfs.append(
            (
                bf.M5DiffBasisFunction(
                    bandname=bandname, nside=nside, fiducial_FWHMEff=fiducial_fwhm
                ),
                m5_weight / 2.0,
            )
        )
        bfs.append(
            (
                bf.M5DiffBasisFunction(
                    bandname=bandname2, nside=nside, fiducial_FWHMEff=fiducial_fwhm
                ),
                m5_weight / 2.0,
            )
        )

    else:
        bfs.append(
            (
                bf.M5DiffBasisFunction(
                    bandname=bandname, nside=nside, fiducial_FWHMEff=fiducial_fwhm
                ),
                m5_weight,
            )
        )

    if bandname2 is not None:
        bfs.append(
            (
                bf.FootprintBasisFunction(
                    bandname=bandname,
                    footprint=footprints,
                    out_of_bounds_val=np.nan,
                    nside=nside,
                ),
                footprint_weight / 2.0,
            )
        )
        bfs.append(
            (
                bf.FootprintBasisFunction(
                    bandname=bandname2,
                    footprint=footprints,
                    out_of_bounds_val=np.nan,
                    nside=nside,
                ),
                footprint_weight / 2.0,
            )
        )
    else:
        bfs.append(
            (
                bf.FootprintBasisFunction(
                    bandname=bandname,
                    footprint=footprints,
                    out_of_bounds_val=np.nan,
                    nside=nside,
                ),
                footprint_weight,
            )
        )

    bfs.append(
        (
            bf.SlewtimeBasisFunction(bandname=bandname, nside=nside),
            slewtime_weight,
        )
    )
    if strict:
        bfs.append((bf.StrictBandBasisFunction(bandname=bandname), stayband_weight))
    else:
        bfs.append((bf.BandChangeBasisFunction(bandname=bandname), stayband_weight))

    if n_obs_template is not None:
        if bandname2 is not None:
            bfs.append(
                (
                    bf.NObsPerYearBasisFunction(
                        bandname=bandname,
                        nside=nside,
                        footprint=footprints.get_footprint(bandname),
                        n_obs=n_obs_template[bandname],
                        season=season,
                        season_start_hour=season_start_hour,
                        season_end_hour=season_end_hour,
                    ),
                    template_weights[bandname] / 2.0,
                )
            )
            bfs.append(
                (
                    bf.NObsPerYearBasisFunction(
                        bandname=bandname2,
                        nside=nside,
                        footprint=footprints.get_footprint(bandname2),
                        n_obs=n_obs_template[bandname2],
                        season=season,
                        season_start_hour=season_start_hour,
                        season_end_hour=season_end_hour,
                    ),
                    template_weights[bandname2] / 2.0,
                )
            )
        else:
            bfs.append(
                (
                    bf.NObsPerYearBasisFunction(
                        bandname=bandname,
                        nside=nside,
                        footprint=footprints.get_footprint(bandname),
                        n_obs=n_obs_template[bandname],
                        season=season,
                        season_start_hour=season_start_hour,
                        season_end_hour=season_end_hour,
                    ),
                    template_weights[bandname],
                )
            )

    # The shared masks
    bandnames = [fn for fn in [bandname, bandname2] if fn is not None]
    bfs.append((bf.BandLoadedBasisFunction(bandnames=bandnames), 0))

    return bfs
