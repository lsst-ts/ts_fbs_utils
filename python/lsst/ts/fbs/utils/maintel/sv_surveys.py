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
import numpy.typing as npt
import rubin_scheduler.scheduler.basis_functions as bf
import rubin_scheduler.scheduler.detailers as detailers
from .ddf_presched import generate_ddf_scheduled_obs
from rubin_scheduler.scheduler.example import simple_greedy_survey
from rubin_scheduler.scheduler.surveys import (
    BlobSurvey,
    GreedySurvey,
    LongGapSurvey,
    ScriptedSurvey,
)
from rubin_scheduler.scheduler.utils import Footprints
from rubin_scheduler.utils import DEFAULT_NSIDE, SURVEY_START_MJD

__all__ = [
    "safety_masks",
    "standard_bf",
    "gen_lvk_templates",
    "gen_template_surveys",
    "blob_for_long",
    "gen_long_gaps_survey",
    "gen_greedy_surveys",
    "generate_blobs",
    "generate_twi_blobs",
    "gen_ddf_surveys",
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
                bf.M5DiffBasisFunction(bandname=bandname, nside=nside),
                m5_weight / 2.0,
            )
        )
        bfs.append(
            (
                bf.M5DiffBasisFunction(bandname=bandname2, nside=nside),
                m5_weight / 2.0,
            )
        )

    else:
        bfs.append((bf.M5DiffBasisFunction(bandname=bandname, nside=nside), m5_weight))

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


def gen_lvk_templates(
    footprints_hp: npt.NDArray,
    nside: int = DEFAULT_NSIDE,
    bands: tuple[str] = ("g", "i"),
    survey_start: float = SURVEY_START_MJD,
    camera_rot_limits: tuple[float, float] = CAMERA_ROT_LIMITS,
    exptime: float = EXPTIME,
    nexp: int = NEXP,
    science_program: str = SCIENCE_PROGRAM,
) -> list[GreedySurvey]:
    """Generate simple area covering single-visit surveys for the end of
    the night, to start to build up early templates in case of LVK alert.

    Parameters
    ----------
    footprints_hp : `np.ndarray` (N, M)
        The (full nside healpix) array containing the footprint values
        for each band used in bands.
    nside : `int`
        The nside for the surveys.
    bands : `list` [ `str` ]
        The bands to run the survey in. Should be g/i, might drop to i.
    survey_start : `float`
        Survey start, in MJD, for the Footprint definition.
    camera_rot_limits : `list` [ `float`, `float` ]
        The camera rotator limits for dithering, in degrees.
    science_program : `str`
        The science_program value to use for the surveys.
    """
    lvk_templates = []
    for band in bands:
        extra_templates_masks = safety_masks(nside=nside, shadow_minutes=15)
        extra_templates_masks.append(
            bf.MaskAfterNObsBasisFunction(nside=nside, n_max=4, bandname=band)
        )
        s = simple_greedy_survey(
            nside=nside,
            bandname=band,
            mask_basis_functions=extra_templates_masks,
            survey_start=survey_start,
            footprints_hp=footprints_hp,
            camera_rot_limits=camera_rot_limits,
            exptime=exptime,
            nexp=nexp,
            science_program=science_program,
            observation_reason=f"template_area_singles_{band}",
        )
        lvk_templates.append(s)

    # Remove detailer (drop this when updated at rubin_scheduler)
    bad_detailer = detailers.Rottep2RotspDesiredDetailer
    for survey in lvk_templates:
        good_dets = []
        for det in survey.detailers:
            if not isinstance(det, bad_detailer):
                good_dets.append(det)
        survey.detailers = good_dets

    return lvk_templates


def gen_template_surveys(
    footprints: Footprints,
    nside: int = DEFAULT_NSIDE,
    band1s: list[str] = ["u", "g", "r", "i", "z", "y"],
    dark_only: list[str] = ["u", "g"],
    ignore_obs: str | list[str] = ["DD", "twilight_near_sun"],
    camera_rot_limits: tuple[float, float] = CAMERA_ROT_LIMITS,
    exptime: float = EXPTIME,
    nexp: int = NEXP,
    u_exptime: float = U_EXPTIME,
    u_nexp: int = U_NEXP,
    n_obs_template: dict | None = None,
    pair_time: float = 33.0,
    area_required: float = 50.0,
    HA_min: float = 2.5,
    HA_max: float = 24 - 2.5,
    max_alt: float = 76.0,
    science_program: str = SCIENCE_PROGRAM,
    blob_survey_params: dict | None = None,
) -> list[BlobSurvey]:
    """Coherent area surveys (BlobSurvey with single visit) that
    are intended to acquire template visits in a convenient yet
    aggressive manner.

    Parameters
    ----------
    footprints : `rubin_scheduler.scheduler.utils.Footprints`
        The Footprints object for the Surveys.
    nside : `int`
        Nside for the surveys.
    band1s : `list` [ `str` ]
        The bands in which to obtain templates, within the Footprints.
    dark_only : `list` [ `str` ]
        The bands to only attempt during dark-time.
    ignore_obs : `str` or `list` [ `str` ]
        Strings to match within scheduler_note to flag observations to ignore.
    camera_rot_limits : `list` [ `float`, `float` ]
        Camera rotator limits (in degrees) for the dither rotation detailer.
    exptime : `float`
        The exposure time for grizy visits.
    nexp : `int`
        The number of exposures per visit for grizy visits.
    u_exptime : `float`
        The exposure time for u band visits.
    u_nexp : `int`
        The number of exposures per visit for u band visits.
    n_obs_template : `dict` { `str` : `int` }
        Number of visits per bandpass before masking the Survey healpix.
        When this many visits are acquired, the pixel is considered "done".
    pair_time : `float`
        The time until the end of the first pass of the Blob.
        Since there is no second filter, this is the amount of time
        spent in the Blob.
    area_required : `float`
        The area required that needs templates, before the BlobSurvey will
        activate.
    HA_min : `float`
        The minimum HA to consider when considering template area.
    HA_max : `float`
        The maximum HA to consider when considering template area.
    max_alt : `float`
        The maximum altitude to use for the Surveys.
        Typically for BlobSurveys this is set lower than the max available,
        to about 76 degrees, to avoid long dome slews near azimuth.
        This is masked separately from the `safety_masks`.
    science_program : `str`
        The science_program to use for visits from these surveys.
    blob_survey_params : `dict` or None
        A dictionary of additional kwargs to pass to the BlobSurvey.
        In particular, the times for typical slews, readtime, etc. are
        useful for setting the number of pointings to schedule within
        pair_time.
    """

    if n_obs_template is None:
        n_obs_template = {"u": 4, "g": 4, "r": 4, "i": 4, "z": 4, "y": 4}

    if blob_survey_params is None:
        blob_survey_params = {
            "slew_approx": 7.5,
            "band_change_approx": 140.0,
            "read_approx": 2.4,
            "flush_time": 30.0,
            "smoothing_kernel": None,
            "nside": nside,
            "seed": 42,
            "dither": "night",
            "twilight_scale": True,
        }

    surveys = []

    for bandname in band1s:
        # Set up Detailers for camera rotator and ordering in altitude.
        detailer_list = []
        detailer_list.append(
            detailers.CameraRotDetailer(
                min_rot=np.min(camera_rot_limits), max_rot=np.max(camera_rot_limits)
            )
        )
        detailer_list.append(detailers.CloseAltDetailer())
        detailer_list.append(
            detailers.BandNexp(bandname="u", nexp=u_nexp, exptime=u_exptime)
        )

        # List to hold tuples of (basis_function_object, weight)
        bfs = []

        bfs.extend(
            standard_bf(
                nside,
                bandname=bandname,
                bandname2=None,
                footprints=footprints,
                n_obs_template=n_obs_template,
            )
        )

        bfs.append(
            (
                bf.AltAzShadowMaskBasisFunction(
                    nside=nside, shadow_minutes=pair_time, max_alt=max_alt, pad=3.0
                ),
                0,
            )
        )

        # not in twilight bf
        bfs.append((bf.TimeToTwilightBasisFunction(time_needed=pair_time), 0.0))
        bfs.append((bf.NotTwilightBasisFunction(), 0.0))
        bfs.append(
            (bf.RevHaMaskBasisFunction(ha_min=HA_min, ha_max=HA_max, nside=nside), 0.0)
        )

        # Need a once in night mask
        bfs.append((bf.NInNightMaskBasisFunction(n_limit=1, nside=nside), 0.0))

        # If u or g, only when moon is down
        if bandname in dark_only:
            bfs.append((bf.MoonAltLimitBasisFunction(alt_limit=-5), 0.0))

        # limit to first year -
        # BUT for SV surveys we know < 1 year, plus "night" is wrong
        # bfs.append((bf.OnlyBeforeNightBasisFunction(night_max=366), 0.0))

        # Mask anything observed n_obs_template times
        bfs.append(
            (
                bf.MaskAfterNObsBasisFunction(
                    nside=nside, n_max=n_obs_template[bandname], bandname=bandname
                ),
                0.0,
            )
        )

        # Add safety masks
        masks = safety_masks(nside=nside, shadow_minutes=pair_time)
        for m in masks:
            bfs.append((m, 0))

        # unpack the basis functions and weights
        weights = [val[1] for val in bfs]
        basis_functions = [val[0] for val in bfs]

        survey_name = "templates, %s" % bandname

        surveys.append(
            BlobSurvey(
                basis_functions,
                weights,
                bandname1=bandname,
                bandname2=None,
                exptime=exptime,
                ideal_pair_time=pair_time,
                survey_name=survey_name,
                science_program=science_program,
                observation_reason=f"template_blob_{bandname}_{pair_time :.1f}",
                ignore_obs=ignore_obs,
                nexp=nexp,
                detailers=detailer_list,
                area_required=area_required,
                **blob_survey_params,
            )
        )

    return surveys


def blob_for_long(
    footprints: Footprints,
    nside: int = DEFAULT_NSIDE,
    band1s: list[str] = ["g"],
    band2s: list[str] = ["i"],
    ignore_obs: str | list[str] = ["DD", "twilight_near_sun", "templates"],
    camera_rot_limits: tuple[float, float] = CAMERA_ROT_LIMITS,
    exptime: float = EXPTIME,
    nexp: int = NEXP,
    u_exptime: float = U_EXPTIME,
    u_nexp: int = U_NEXP,
    n_obs_template: dict = None,
    pair_time: float = 33.0,
    season: float = 365.25,
    season_start_hour: float = -4.0,
    season_end_hour: float = 2.0,
    HA_min: float = 12,
    HA_max: float = 24 - 3.5,
    max_alt: float = 76.0,
    m5_weight: float = 6.0,
    footprint_weight: float = 1.5,
    slewtime_weight: float = 3.0,
    stayband_weight: float = 3.0,
    template_weight: float = 12.0,
    u_template_weight: float = 50.0,
    g_template_weight: float = 50.0,
    night_pattern: list[bool] = [True, True],
    time_after_twi: float = 30.0,
    blob_names: list[str] = [],
    scheduled_respect: float = 30.0,
    science_program: str = SCIENCE_PROGRAM,
    observation_reason: str | None = None,
    blob_survey_params: dict | None = None,
) -> list[BlobSurvey]:
    """
    Generate surveys that take observations in blobs, specifically for the
    long-gaps (triplets) survey.

    Parameters
    ----------
    footprints : `rubin_scheduler.scheduler.utils.Footprints`
        The Footprints object for the Surveys.
    nside : `int`
        The HEALpix nside to use. Default to DEFAULT_NSIDE.
    band1s : `list` [`str`]
        The bandnames for the first band in a pair.
        Default ["g"].
    band2s : `list` of `str`
        The band names for the second in the pair (None if unpaired).
        Default ["i"].
    ignore_obs : `str` or `list` of `str`
        Ignore observations by surveys that include the given substring(s).
    camera_rot_limits : `list` of `float`
        The limits to impose when rotationally dithering the camera (degrees).
        Default [-80., 80.].
    exptime : `float`
        The exposure time for grizy visits.
    nexp : `int`
        The number of exposures per visit for grizy visits.
    u_exptime : `float`
        The exposure time for u band visits.
    u_nexp : `int`
        The number of exposures per visit for u band visits.
    n_obs_template : `dict`
        The number of observations to take every season in each band.
        If None, sets to 3 each. Default None.
    pair_time : `float`
        The ideal time between pairs (minutes). Default 33.
    season : float
        The length of season (i.e., how long before templates expire) (days)
    season_start_hour : `float`
        Hour angle limits to use when gathering templates.
        Default -4 (hours)
    season_end_hour : `float`
       Hour angle limits to use when gathering templates.
       Default +2 (hours)
    HA_min : `float`
        The minimum HA to consider when considering template area.
    HA_max : `float`
        The maximum HA to consider when considering template area.
    max_alt : `float`
        The maximum altitude to use for the Surveys.
        Typically for BlobSurveys this is set lower than the max available,
        to about 76 degrees, to avoid long dome slews near azimuth.
        This is masked separately from the `safety_masks`.
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
        are so few u-visits, it can be helpful to turn this up a
        little higher than the standard template_weight kwarg.
    g_template_weight : `float`
        The weight to place on getting image templates in g-band.
    night_pattern : `list` [ `bool` ]
        Flips the blob_for_long surveys on or off, on a pattern of
        nights. These surveys typically don't execute every night.
    time_after_twi : `float`
        Don't start before this many minutes pass after -18 degree twilight.
    blob_names : `list` [ `str` ]
        Strings to match in scheduler_note to ensure surveys don't execute
        more than once per night.
    scheduled_respect : `float`
        Ensure that blobs don't start within this many minutes of scheduled
        observations (from a ScriptedSurvey).
    science_program : `str`
        The science_program to use for visits from these surveys.
    observation_reason : `str` or None
        The value to use for the observation_reason.
    blob_survey_params : `dict` or None
        A dictionary of additional kwargs to pass to the BlobSurvey.
        In particular, the times for typical slews, readtime, etc. are
        useful for setting the number of pointings to schedule within
        pair_time.
    """

    if blob_survey_params is None:
        blob_survey_params = {
            "slew_approx": 7.5,
            "band_change_approx": 140.0,
            "read_approx": 2.4,
            "flush_time": 30.0,
            "smoothing_kernel": None,
            "nside": nside,
            "seed": 42,
            "dither": "night",
            "twilight_scale": True,
        }

    surveys = []
    if n_obs_template is None:
        n_obs_template = {"u": 3, "g": 3, "r": 3, "i": 3, "z": 3, "y": 3}

    times_needed = [pair_time, pair_time * 2]
    for bandname, bandname2 in zip(band1s, band2s):
        detailer_list = []
        detailer_list.append(
            detailers.CameraRotDetailer(
                min_rot=np.min(camera_rot_limits), max_rot=np.max(camera_rot_limits)
            )
        )
        detailer_list.append(detailers.CloseAltDetailer())
        if "u" in [bandname, bandname2]:
            detailer_list.append(
                detailers.BandNexp(bandname="u", nexp=u_nexp, exptime=u_exptime)
            )

        # List to hold tuples of (basis_function_object, weight)
        bfs = []

        bfs.extend(
            standard_bf(
                nside,
                bandname=bandname,
                bandname2=bandname2,
                m5_weight=m5_weight,
                footprint_weight=footprint_weight,
                slewtime_weight=slewtime_weight,
                stayband_weight=stayband_weight,
                template_weight=template_weight,
                u_template_weight=u_template_weight,
                g_template_weight=g_template_weight,
                footprints=footprints,
                n_obs_template=n_obs_template,
                season=season,
                season_start_hour=season_start_hour,
                season_end_hour=season_end_hour,
            )
        )

        # Make sure we respect scheduled observations
        bfs.append((bf.TimeToScheduledBasisFunction(time_needed=scheduled_respect), 0))

        # Masks, give these 0 weight
        shadow_minutes = pair_time * 2 + 5
        bfs.append(
            (
                bf.AltAzShadowMaskBasisFunction(
                    nside=nside, shadow_minutes=shadow_minutes, max_alt=max_alt, pad=3.0
                ),
                0.0,
            )
        )
        if bandname2 is None:
            time_needed = times_needed[0]
        else:
            time_needed = times_needed[1]
        bfs.append((bf.TimeToTwilightBasisFunction(time_needed=time_needed), 0.0))
        bfs.append((bf.NotTwilightBasisFunction(), 0.0))
        bfs.append((bf.AfterEveningTwiBasisFunction(time_after=time_after_twi), 0.0))
        bfs.append(
            (bf.HaMaskBasisFunction(ha_min=HA_min, ha_max=HA_max, nside=nside), 0.0)
        )
        # don't execute every night
        bfs.append((bf.NightModuloBasisFunction(night_pattern), 0.0))
        # only execute one blob per night
        bfs.append((bf.OnceInNightBasisFunction(notes=blob_names), 0))

        # Add safety masks
        masks = safety_masks(nside=nside, shadow_minutes=shadow_minutes)
        for m in masks:
            bfs.append((m, 0))

        # unpack the basis functions and weights
        weights = [val[1] for val in bfs]
        basis_functions = [val[0] for val in bfs]
        if bandname2 is None:
            survey_name = "blob_long, %s" % bandname
        else:
            survey_name = "blob_long, %s%s" % (bandname, bandname2)
        if bandname2 is not None:
            detailer_list.append(detailers.TakeAsPairsDetailer(bandname=bandname2))

        if observation_reason is None:
            observation_reason = f"triplet_pairs_{bandname}{bandname2}_{pair_time :.1f}"

        surveys.append(
            BlobSurvey(
                basis_functions,
                weights,
                bandname1=bandname,
                bandname2=bandname2,
                exptime=exptime,
                ideal_pair_time=pair_time,
                survey_name=survey_name,
                ignore_obs=ignore_obs,
                nexp=nexp,
                detailers=detailer_list,
                science_program=science_program,
                observation_reason=observation_reason,
                **blob_survey_params,
            )
        )

    return surveys


def gen_long_gaps_survey(
    footprints: Footprints,
    nside: int = DEFAULT_NSIDE,
    camera_rot_limits: tuple[float, float] = CAMERA_ROT_LIMITS,
    exptime: float = EXPTIME,
    nexp: int = NEXP,
    u_exptime: float = U_EXPTIME,
    u_nexp: int = U_NEXP,
    pair_time: float = 33.0,
    night_pattern: list[bool] = [True, True],
    gap_range: list[float] = [2, 7],
    HA_min: float = 12,
    HA_max: float = 24 - 3.5,
    time_after_twi: float = 120,
    u_template_weight: float = 50.0,
    g_template_weight: float = 50.0,
    science_program: str = SCIENCE_PROGRAM,
    blob_survey_params: dict | None = None,
) -> list[LongGapSurvey]:
    """
    Generate long-gaps (triplets) surveys.

    Parameters
    -----------
    footprints : `rubin_scheduler.scheduler.utils.footprints.Footprints`
        The footprints to be used for the long-gaps surveys.
    nside : `int`
        The nside for the surveys.
    camera_rot_limits : `list` of `float`
        The limits to impose when rotationally dithering the camera (degrees).
        Default [-80., 80.].
    exptime : `float`
        The exposure time for grizy visits.
    nexp : `int`
        The number of exposures per visit for grizy visits.
    u_exptime : `float`
        The exposure time for u band visits.
    u_nexp : `int`
        The number of exposures per visit for u band visits.
    pair_time : `float`
        The ideal time between pairs (minutes). Default 33.
    night_pattern : `list` [ `bool` ]
        Which nights to let the survey execute.
        Default of [True, True] executes every night.
    gap_range : `list` [ `float` ]
        Range of times to attempt to gather pairs (hours).
        Default [2, 7].
    HA_min : `float`
        The hour angle limits passed to the initial blob scheduler.
        Default 12 (hours)
    HA_max : `float`
        The hour angle limits passed to the initial blob scheduler.
        Default 20.5 (hours).
    time_after_twi : `float`
        The time after evening twilight to attempt long gaps (minutes).
        Default 120.
    u_template_weight : `float`
        The weight to place on getting image templates in u-band. Since there
        are so few u-visits, it can be helpful to turn this up a
        little higher than the standard template_weight kwarg.
    g_template_weight : `float`
        The weight to place on getting image templates in u-band. Since there
        are so few u-visits, it can be helpful to turn this up a
        little higher than the standard template_weight kwarg.
    science_program : `str`
        The science_program to use for visits from these surveys.
    blob_survey_params : `dict` or None
        A dictionary of additional kwargs to pass to the BlobSurvey.
        In particular, the times for typical slews, readtime, etc. are
        useful for setting the number of pointings to schedule within
        pair_time.
    """

    surveys = []
    f1 = ["g", "r", "i"]
    f2 = ["r", "i", "z"]
    # Maybe force scripted to not go in twilight?
    blob_names = []
    for fn1, fn2 in zip(f1, f2):
        for ab in ["a", "b"]:
            blob_names.append("blob_long, %s%s, %s" % (fn1, fn2, ab))
    for bandname1, bandname2 in zip(f1, f2):
        blob = blob_for_long(
            footprints=footprints,
            camera_rot_limits=camera_rot_limits,
            exptime=exptime,
            nexp=nexp,
            u_exptime=u_exptime,
            u_nexp=u_nexp,
            pair_time=pair_time,
            nside=nside,
            band1s=[bandname1],
            band2s=[bandname2],
            night_pattern=night_pattern,
            time_after_twi=time_after_twi,
            HA_min=HA_min,
            HA_max=HA_max,
            u_template_weight=u_template_weight,
            g_template_weight=g_template_weight,
            blob_names=blob_names,
            science_program=science_program,
            blob_survey_params=blob_survey_params,
        )
        # If one of the bands is u, need to add detailer
        scripted = ScriptedSurvey(
            safety_masks(nside, shadow_minutes=pair_time),
            nside=nside,
            ignore_obs=["blob", "DDF", "twi", "pair", "templates"],
            science_program=science_program,
        )
        surveys.append(
            LongGapSurvey(blob[0], scripted, gap_range=gap_range, avoid_zenith=True)
        )

    return surveys


def gen_greedy_surveys(
    nside: int = DEFAULT_NSIDE,
    bands: list[str] = ["r", "i", "z", "y"],
    ignore_obs: list[str] = ["DD", "twilight_near_sun, templates"],
    camera_rot_limits: tuple[float, float] = CAMERA_ROT_LIMITS,
    exptime: float = EXPTIME,
    nexp: int = NEXP,
    u_exptime: float = U_EXPTIME,
    u_nexp: int = U_NEXP,
    shadow_minutes: float = 15.0,
    max_alt: float = 76.0,
    m5_weight: float = 3.0,
    footprint_weight: float = 0.75,
    slewtime_weight: float = 3.0,
    stayband_weight: float = 100.0,
    repeat_weight: float = -1.0,
    footprints: Footprints | None = None,
    science_program=SCIENCE_PROGRAM,
) -> list[GreedySurvey]:
    """Generate greedy (single-best choice visits) Surveys.

    Parameters
    ----------
    nside : `int`
        The HEALpix nside to use
    bands : `list` [ `str` ]
        Bands in which to generate greedy surveys.
        Default ['r', 'i', 'z', 'y'].
    ignore_obs : `str` or `list` of `str`
        Ignore observations by surveys that include the given substring(s).
    camera_rot_limits : `list` [ `float` ]
        The limits to impose when rotationally dithering the camera (degrees).
        Default [-80., 80.].
    exptime : `float`
        The exposure time for grizy visits.
    nexp : `int`
        The number of exposures per visit for grizy visits.
    u_exptime : `float`
        The exposure time for u band visits.
    u_nexp : `int`
        The number of exposures per visit for u band visits.
    shadow_minutes : `float`
        Used to mask regions around zenith (minutes).
    max_alt : `float`
        The maximium altitude to use when masking zenith (degrees).
    m5_weight : `float`
        The weight for the 5-sigma depth difference basis function.
    footprint_weight : `float`
        The weight on the survey footprint basis function.
    slewtime_weight : `float`
        The weight on the slewtime basis function.
    stayband_weight : `float`
        The weight on basis function that tries to stay avoid band changes.
    repeat_weight : `float`
        Weight that enhances (+ values) or decreases (- values) the likelihood
        of revisiting the same pointing within a two-hour time gap.
    """
    # Define the extra parameters that are used in the greedy survey.
    greed_survey_params = {
        "block_size": 1,
        "smoothing_kernel": None,
        "seed": 42,
        "camera": "LSST",
        "dither": "night",
    }

    surveys = []
    detailer_list = [
        detailers.CameraRotDetailer(
            min_rot=np.min(camera_rot_limits), max_rot=np.max(camera_rot_limits)
        )
    ]
    if "u" in bands:
        detailer_list.append(
            detailers.BandNexp(bandname="u", nexp=u_nexp, exptime=u_exptime)
        )

    for bandname in bands:
        bfs = []
        bfs.extend(
            standard_bf(
                nside,
                bandname=bandname,
                bandname2=None,
                m5_weight=m5_weight,
                footprint_weight=footprint_weight,
                slewtime_weight=slewtime_weight,
                stayband_weight=stayband_weight,
                template_weight=0,
                u_template_weight=0,
                g_template_weight=0,
                footprints=footprints,
                n_obs_template=None,
                strict=False,
            )
        )

        bfs.append(
            (
                bf.VisitRepeatBasisFunction(
                    gap_min=0, gap_max=2 * 60.0, bandname=None, nside=nside, npairs=20
                ),
                repeat_weight,
            )
        )
        # Masks, give these 0 weight
        bfs.append(
            (
                bf.AltAzShadowMaskBasisFunction(
                    nside=nside, shadow_minutes=shadow_minutes, max_alt=max_alt, pad=3.0
                ),
                0,
            )
        )

        masks = safety_masks(nside=nside, shadow_minutes=shadow_minutes)
        for m in masks:
            bfs.append((m, 0))

        weights = [val[1] for val in bfs]
        basis_functions = [val[0] for val in bfs]
        surveys.append(
            GreedySurvey(
                basis_functions,
                weights,
                exptime=exptime,
                bandname=bandname,
                nside=nside,
                ignore_obs=ignore_obs,
                nexp=nexp,
                detailers=detailer_list,
                survey_name=f"greedy {bandname}",
                science_program=science_program,
                observation_reason=f"singles_{bandname}",
                **greed_survey_params,
            )
        )

    return surveys


def generate_blobs(
    footprints: Footprints,
    nside: int = DEFAULT_NSIDE,
    band1s: list[str] = ["u", "u", "g", "r", "i", "z", "y"],
    band2s: list[str] = ["g", "r", "r", "i", "z", "y", "y"],
    ignore_obs: str | list[str] = ["DD", "twilight_near_sun", "templates"],
    camera_rot_limits: tuple[float, float] = CAMERA_ROT_LIMITS,
    exptime: float = EXPTIME,
    nexp: int = NEXP,
    u_exptime: float = U_EXPTIME,
    u_nexp: int = U_NEXP,
    n_obs_template: dict = None,
    pair_time: float = 33.0,
    season: float = 365.25,
    season_start_hour: float = -4.0,
    season_end_hour: float = 2.0,
    max_alt: float = 76.0,
    m5_weight: float = 6.0,
    footprint_weight: float = 1.5,
    slewtime_weight: float = 3.0,
    stayband_weight: float = 3.0,
    template_weight: float = 12.0,
    u_template_weight: float = 50.0,
    g_template_weight: float = 50.0,
    repeat_weight: float = -20,
    good_seeing: dict = {"g": 3, "r": 3, "i": 3},
    good_seeing_weight: float = 3.0,
    survey_start: float = SURVEY_START_MJD,
    scheduled_respect: float = 45.0,
    science_program: str = SCIENCE_PROGRAM,
    blob_survey_params: dict | None = None,
) -> list[BlobSurvey]:
    """
    Generate surveys that take observations in blobs.

    Parameters
    ----------
    footprints : `rubin_scheduler.scheduler.utils.Footprints`
        The Footprints object for the Surveys.
    nside : `int`
        The HEALpix nside to use. Default to DEFAULT_NSIDE.
    band1s : `list` [`str`]
        The bandnames for the first band in a pair.
    band2s : `list` of `str`
        The band names for the second in the pair (None if unpaired).
    ignore_obs : `str` or `list` of `str`
        Ignore observations by surveys that include the given substring(s).
    camera_rot_limits : `list` of `float`
        The limits to impose when rotationally dithering the camera (degrees).
    exptime : `float`
        The exposure time for grizy visits.
    nexp : `int`
        The number of exposures per visit for grizy visits.
    u_exptime : `float`
        The exposure time for u band visits.
    u_nexp : `int`
        The number of exposures per visit for u band visits.
    n_obs_template : `dict`
        The number of observations to take every season in each band.
        If None, sets to 3 each. Default None.
    pair_time : `float`
        The ideal time between pairs (minutes). Default 33.
    season : float
        The length of season (i.e., how long before templates expire) (days)
    season_start_hour : `float`
        Hour angle limits to use when gathering templates.
        Default -4 (hours)
    season_end_hour : `float`
       Hour angle limits to use when gathering templates.
       Default +2 (hours)
    max_alt : `float`
        The maximum altitude to use for the Surveys.
        Typically for BlobSurveys this is set lower than the max available,
        to about 76 degrees, to avoid long dome slews near azimuth.
        This is masked separately from the `safety_masks`.
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
        are so few u-visits, it can be helpful to turn this up a
        little higher than the standard template_weight kwarg.
    g_template_weight : `float`
        The weight to place on getting image templates in g-band.
    repeat_weight : `float`
        Weight that enhances (+ values) or decreases (- values) the likelihood
        of revisiting the same pointing within a two-hour time gap.
    good_seeing : `dict`
        Number of good-seeing (<6") images per band for basis function
        directing good-atmospheric seeing into particular bands.
    good_seeing_weight : `float`
        The weight to place on getting good-seeing images in the
        good_seeing bands.
    survey_start : `float`
        The mjd that the survey started (used for determining season for
        counting good seeing images within a season).
    scheduled_respect : `float`
        Ensure that blobs don't start within this many minutes of scheduled
        observations (from a ScriptedSurvey).
    science_program : `str`
        The science_program to use for visits from these surveys.
    blob_survey_params : `dict` or None
        A dictionary of additional kwargs to pass to the BlobSurvey.
        In particular, the times for typical slews, readtime, etc. are
        useful for setting the number of pointings to schedule within
        pair_time.
    """

    if blob_survey_params is None:
        blob_survey_params = {
            "slew_approx": 7.5,
            "band_change_approx": 140.0,
            "read_approx": 2.4,
            "flush_time": 30.0,
            "smoothing_kernel": None,
            "nside": nside,
            "seed": 42,
            "dither": "night",
            "twilight_scale": True,
        }

    if n_obs_template is None:
        n_obs_template = {"u": 3, "g": 3, "r": 3, "i": 3, "z": 3, "y": 3}

    surveys = []

    times_needed = [pair_time, pair_time * 2]
    for bandname, bandname2 in zip(band1s, band2s):
        detailer_list = []
        detailer_list.append(
            detailers.CameraRotDetailer(
                min_rot=np.min(camera_rot_limits), max_rot=np.max(camera_rot_limits)
            )
        )
        detailer_list.append(detailers.CloseAltDetailer())
        detailer_list.append(detailers.FlushForSchedDetailer())
        if bandname == "u" or bandname2 == "u":
            detailer_list.append(
                detailers.BandNexp(bandname="u", nexp=u_nexp, exptime=u_exptime)
            )

        # List to hold tuples of (basis_function_object, weight)
        bfs = []

        bfs.extend(
            standard_bf(
                nside,
                bandname=bandname,
                bandname2=bandname2,
                m5_weight=m5_weight,
                footprint_weight=footprint_weight,
                slewtime_weight=slewtime_weight,
                stayband_weight=stayband_weight,
                template_weight=template_weight,
                u_template_weight=u_template_weight,
                g_template_weight=g_template_weight,
                footprints=footprints,
                n_obs_template=n_obs_template,
                season=season,
                season_start_hour=season_start_hour,
                season_end_hour=season_end_hour,
            )
        )

        bfs.append(
            (
                bf.VisitRepeatBasisFunction(
                    gap_min=0, gap_max=3 * 60.0, bandname=None, nside=nside, npairs=20
                ),
                repeat_weight,
            )
        )

        # Insert things for getting good seeing templates
        if bandname2 is not None:
            if bandname in list(good_seeing.keys()):
                bfs.append(
                    (
                        bf.NGoodSeeingBasisFunction(
                            bandname=bandname,
                            nside=nside,
                            mjd_start=survey_start,
                            footprint=footprints.get_footprint(bandname),
                            n_obs_desired=good_seeing[bandname],
                        ),
                        good_seeing_weight,
                    )
                )
            if bandname2 in list(good_seeing.keys()):
                bfs.append(
                    (
                        bf.NGoodSeeingBasisFunction(
                            bandname=bandname2,
                            nside=nside,
                            mjd_start=survey_start,
                            footprint=footprints.get_footprint(bandname2),
                            n_obs_desired=good_seeing[bandname2],
                        ),
                        good_seeing_weight,
                    )
                )
        else:
            if bandname in list(good_seeing.keys()):
                bfs.append(
                    (
                        bf.NGoodSeeingBasisFunction(
                            bandname=bandname,
                            nside=nside,
                            mjd_start=survey_start,
                            footprint=footprints.get_footprint(bandname),
                            n_obs_desired=good_seeing[bandname],
                        ),
                        good_seeing_weight,
                    )
                )
        # Make sure we respect scheduled observations
        bfs.append((bf.TimeToScheduledBasisFunction(time_needed=scheduled_respect), 0))
        # Masks, give these 0 weight

        shadow_minutes = pair_time * 2 + 5
        bfs.append(
            (
                bf.AltAzShadowMaskBasisFunction(
                    nside=nside, shadow_minutes=shadow_minutes, max_alt=max_alt, pad=3.0
                ),
                0.0,
            )
        )
        if bandname2 is None:
            time_needed = times_needed[0]
        else:
            time_needed = times_needed[1]
        bfs.append((bf.TimeToTwilightBasisFunction(time_needed=time_needed), 0.0))
        bfs.append((bf.NotTwilightBasisFunction(), 0.0))

        # Add safety masks
        masks = safety_masks(nside=nside, shadow_minutes=shadow_minutes)
        for m in masks:
            bfs.append((m, 0))

        # unpack the basis functions and weights
        weights = [val[1] for val in bfs]
        basis_functions = [val[0] for val in bfs]

        # Set survey_name
        if bandname2 is None:
            survey_name = "pair_%i, %s" % (pair_time, bandname)
        else:
            survey_name = "pair_%i, %s%s" % (pair_time, bandname, bandname2)
        if bandname2 is not None:
            detailer_list.append(detailers.TakeAsPairsDetailer(bandname=bandname2))

        observation_reason = f"pairs_{bandname}"
        if bandname2 is not None:
            observation_reason += f"{bandname2}"
        observation_reason += f"_{pair_time :.1f}"

        surveys.append(
            BlobSurvey(
                basis_functions,
                weights,
                bandname1=bandname,
                bandname2=bandname2,
                exptime=exptime,
                ideal_pair_time=pair_time,
                survey_name=survey_name,
                science_program=science_program,
                observation_reason=observation_reason,
                ignore_obs=ignore_obs,
                nexp=nexp,
                detailers=detailer_list,
                **blob_survey_params,
            )
        )

    return surveys


def generate_twi_blobs(
    footprints: Footprints,
    nside: int = DEFAULT_NSIDE,
    band1s: list[str] = ["r", "i", "z", "y"],
    band2s: list[str] = ["i", "z", "y", "y"],
    ignore_obs: str | list[str] = ["DD", "twilight_near_sun", "templates"],
    camera_rot_limits: tuple[float, float] = CAMERA_ROT_LIMITS,
    exptime: float = EXPTIME,
    nexp: int = NEXP,
    n_obs_template: dict = None,
    pair_time: float = 15.0,
    season: float = 365.25,
    season_start_hour: float = -4.0,
    season_end_hour: float = 2.0,
    max_alt: float = 76.0,
    m5_weight: float = 6.0,
    footprint_weight: float = 1.5,
    slewtime_weight: float = 3.0,
    stayband_weight: float = 3.0,
    template_weight: float = 12.0,
    repeat_weight: float = -1,
    scheduled_respect: float = 15.0,
    night_pattern: list[bool] | None = None,
    science_program: str = SCIENCE_PROGRAM,
    blob_survey_params: dict | None = None,
) -> list[BlobSurvey]:
    """
    Generate surveys that take observations in blobs, for twilight time.
    Shorter blobs, different weights for the basis functions.

    Parameters
    ----------
    footprints : `rubin_scheduler.scheduler.utils.Footprints`
        The Footprints object for the Surveys.
    nside : `int`
        The HEALpix nside to use. Default to DEFAULT_NSIDE.
    band1s : `list` [`str`]
        The bandnames for the first band in a pair.
    band2s : `list` [ `str` ]
        The band names for the second in the pair (None if unpaired).
    ignore_obs : `str` or `list` [ `str` ]
        Ignore observations by surveys that include the given substring(s).
    camera_rot_limits : `list` [ `float` ]
        The limits to impose when rotationally dithering the camera (degrees).
    exptime : `float`
        Exposure time for visits.
    nexp : `int`
        Number of exposures per visit.
    n_obs_template : `dict`
        The number of observations to take every season in each band.
        If None, sets to 3 each. Default None.
    pair_time : `float`
        The ideal time between pairs (minutes). Default 33.
    season : float
        The length of season (i.e., how long before templates expire) (days)
    season_start_hour : `float`
        Hour angle limits to use when gathering templates.
        Default -4 (hours)
    season_end_hour : `float`
       Hour angle limits to use when gathering templates.
       Default +2 (hours)
    max_alt : `float`
        The maximum altitude to use for the Surveys.
        Typically for BlobSurveys this is set lower than the max available,
        to about 76 degrees, to avoid long dome slews near azimuth.
        This is masked separately from the `safety_masks`.
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
    repeat_weight : `float`
        Weight that enhances (+ values) or decreases (- values) the likelihood
        of revisiting the same pointing within a two-hour time gap.
    scheduled_respect : `float`
        Ensure that blobs don't start within this many minutes of scheduled
        observations (from a ScriptedSurvey).
    night_pattern : `list` [ `bool` ]
        Which nights to let the survey execute (should be the opposite of
        the pattern for the NEO twilight survey).
        Default of [True, True] executes every night.
    science_program : `str`
        The science_program to use for visits from these surveys.
    blob_survey_params : `dict` or None
        A dictionary of additional kwargs to pass to the BlobSurvey.
        In particular, the times for typical slews, readtime, etc. are
        useful for setting the number of pointings to schedule within
        pair_time.
    """

    if blob_survey_params is None:
        blob_survey_params = {
            "slew_approx": 7.5,
            "band_change_approx": 140.0,
            "read_approx": 2.4,
            "flush_time": 30.0,
            "smoothing_kernel": None,
            "nside": nside,
            "seed": 42,
            "dither": "night",
            "twilight_scale": True,
        }

    surveys = []

    if n_obs_template is None:
        n_obs_template = {"u": 3, "g": 3, "r": 3, "i": 3, "z": 3, "y": 3}

    times_needed = [pair_time, pair_time * 2]
    for bandname, bandname2 in zip(band1s, band2s):
        detailer_list = []
        detailer_list.append(
            detailers.CameraRotDetailer(
                min_rot=np.min(camera_rot_limits), max_rot=np.max(camera_rot_limits)
            )
        )
        detailer_list.append(detailers.CloseAltDetailer())
        detailer_list.append(detailers.FlushForSchedDetailer())
        # List to hold tuples of (basis_function_object, weight)
        bfs = []

        bfs.extend(
            standard_bf(
                nside,
                bandname=bandname,
                bandname2=bandname2,
                m5_weight=m5_weight,
                footprint_weight=footprint_weight,
                slewtime_weight=slewtime_weight,
                stayband_weight=stayband_weight,
                template_weight=template_weight,
                u_template_weight=0,
                g_template_weight=0,
                footprints=footprints,
                n_obs_template=n_obs_template,
                season=season,
                season_start_hour=season_start_hour,
                season_end_hour=season_end_hour,
            )
        )

        bfs.append(
            (
                bf.VisitRepeatBasisFunction(
                    gap_min=0, gap_max=2 * 60.0, bandname=None, nside=nside, npairs=20
                ),
                repeat_weight,
            )
        )

        # Make sure we respect scheduled observations
        bfs.append((bf.TimeToScheduledBasisFunction(time_needed=scheduled_respect), 0))
        # Masks, give these 0 weight
        shadow_minutes = pair_time * 2 + 5
        bfs.append(
            (
                bf.AltAzShadowMaskBasisFunction(
                    nside=nside,
                    shadow_minutes=shadow_minutes,
                    max_alt=max_alt,
                    pad=3.0,
                ),
                0.0,
            )
        )
        if bandname2 is None:
            time_needed = times_needed[0]
        else:
            time_needed = times_needed[1]
        bfs.append(
            (bf.TimeToTwilightBasisFunction(time_needed=time_needed, alt_limit=12), 0.0)
        )

        # Let's turn off twilight blobs on nights where we are
        # doing NEO hunts
        bfs.append((bf.NightModuloBasisFunction(pattern=night_pattern), 0))

        # Add safety masks
        masks = safety_masks(nside=nside, shadow_minutes=shadow_minutes)
        for m in masks:
            bfs.append((m, 0))

        # unpack the basis functions and weights
        weights = [val[1] for val in bfs]
        basis_functions = [val[0] for val in bfs]

        # Set survey names
        if bandname2 is None:
            survey_name = "pair_%i, %s" % (pair_time, bandname)
        else:
            survey_name = "pair_%i, %s%s" % (pair_time, bandname, bandname2)

        if bandname2 is not None:
            detailer_list.append(detailers.TakeAsPairsDetailer(bandname=bandname2))

        observation_reason = f"pairs_{bandname}"
        if bandname2 is not None:
            observation_reason += f"{bandname2}"
        observation_reason = +f"_{pair_time :.1f}"

        surveys.append(
            BlobSurvey(
                basis_functions,
                weights,
                bandname1=bandname,
                bandname2=bandname2,
                exptime=exptime,
                ideal_pair_time=pair_time,
                survey_name=survey_name,
                science_program=science_program,
                observation_reason=observation_reason,
                ignore_obs=ignore_obs,
                nexp=nexp,
                detailers=detailer_list,
                **blob_survey_params,
            )
        )

    return surveys


def gen_ddf_surveys(
    ddf_config_file: str,
    detailer_list: list[detailers.BaseDetailer] | None = None,
    nside: int = DEFAULT_NSIDE,
    expt: dict | None = None,
    nexp: dict | None = None,
    survey_start: float = SURVEY_START_MJD,
    survey_length: int = 10,
    science_program: str = SCIENCE_PROGRAM,
) -> list[ScriptedSurvey]:
    """Generate surveys for DDF observations.

    Parameters
    ----------
    ddf_config_file : `str`
        Filename for the DDF configuration information.
    detailer_list : `list` [ `rubin_scheduler.scheduler.Detailer` ]
        Detailers for DDFs. Default None.
    nside : `int`
        Nside for the survey. Used for mask basis functions.
    expt : `dict`  { `str` : `float` } or None
        Exposure time for DDF visits.
        Default of None uses defaults of EXPTIME/U_EXPTIME.
    nexp : `dict` { `str` : `int` } or None
        Number of exposures per visit.
        Default of None uses defaults of NEXP/U_NEXP.
    survey_start : `float`
        Start MJD of the survey. Used for prescheduling DDF visits.
    survey_length : `float`
        Length of the survey. Used for prescheduling DDF visits.
    science_program : `str`
        Name of the science program for the Survey.
    """
    if expt is None:
        expt = {
            "u": U_EXPTIME,
            "g": EXPTIME,
            "r": EXPTIME,
            "i": EXPTIME,
            "z": EXPTIME,
            "y": EXPTIME,
        }
    if nexp is None:
        nexp = {"u": U_NEXP, "g": NEXP, "r": NEXP, "i": NEXP, "z": NEXP, "y": NEXP}

    if detailer_list is None:
        dither_detailer = detailers.DitherDetailer(per_night=True, max_dither=0.2)

        dither_detailer = detailers.SplitDetailer(
            dither_detailer, detailers.EuclidDitherDetailer()
        )
        detailer_list = [
            detailers.CameraRotDetailer(min_rot=-75, max_rot=75),
            dither_detailer,
            detailers.BandSortDetailer(),
        ]
    # Append a u-band exposure time detailer regardless (override)
    u_detailer = detailers.BandNexp(bandname="u", nexp=nexp["u"], exptime=expt["u"])
    detailer_list.append(u_detailer)

    obs_array = generate_ddf_scheduled_obs(
        ddf_config_file,
        expt=expt,
        nsnaps=nexp,
        mjd_start=survey_start,
        survey_length=survey_length,
        sun_alt_max=-24,
        science_program=science_program,
    )

    survey1 = ScriptedSurvey(
        safety_masks(nside, shadow_minutes=30),
        nside=nside,
        detailers=detailer_list,
        survey_name="deep drilling",
    )
    survey1.set_script(obs_array)

    return [survey1]
