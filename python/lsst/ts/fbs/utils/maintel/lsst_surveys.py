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

__all__ = (
    "standard_masks",
    "standard_bf",
    "gen_long_gaps_survey",
    "gen_template_surveys",
    "gen_greedy_surveys",
    "generate_blobs",
    "generate_twilight_near_sun",
)

import copy
from typing import Any

import numpy as np
import numpy.typing as npt
import rubin_scheduler.scheduler.basis_functions as bf
import rubin_scheduler.scheduler.detailers as detailers
from rubin_scheduler.scheduler.surveys import (
    BlobPairsSurvey,
    BlobSurvey,
    GreedySurvey,
    LongGapSurvey,
    ScriptedSurvey,
)
from rubin_scheduler.scheduler.utils import ConstantFootprint, Footprints, ecliptic_area
from rubin_scheduler.utils import (
    DEFAULT_NSIDE,
    SURVEY_START_MJD,
    declination_dependent_fwhm,
)

# Set up values to use as kwarg defaults.
EXPTIME = 30.0
U_EXPTIME = 38.0
SEEING_FWHM_MAX_ZENITH_DEFAULT = {
    "u": 1.0,
    "g": 1.0,
    "r": 1.0,
    "i": 1.0,
    "z": 1.0,
    "y": 1.0,
}
CAMERA_ROT_LIMITS = (-80.0, 80.0)
SCIENCE_PROGRAM = "BLOCK-430"

BLOB_SURVEY_PARAMS_DEFAULTS = {
    "slew_approx": 8.0,
    "band_change_approx": 200.0,
    "read_approx": 3.07,
    "flush_time": 30.0,
    "smoothing_kernel": None,
    "nside": DEFAULT_NSIDE,
    "seed": 42,
    "twilight_scale": True,
    "check_scheduled": True,
}

STANDARD_MASK_DEFAULTS = {
    "nside": DEFAULT_NSIDE,
    "wind_speed_maximum": 40,
    "min_alt": 20,
    "max_alt": 86.5,
    "shadow_minutes": 2,
    "apply_cloud_mask": True,
    "cloud_limit": 8,
    "apply_time_limited_shadow": False,
    "time_to_sunrise": 3.0,
    "min_az_sunrise": 150,
    "max_az_sunrise": 250,
}


def standard_masks(
    nside: int = DEFAULT_NSIDE,
    moon_distance: float = 30.0,
    wind_speed_maximum: float = 20.0,
    min_alt: float = 20.0,
    max_alt: float = 86.5,
    min_az: float = 0.0,
    max_az: float = 360.0,
    shadow_minutes: float = 0.0,
    apply_time_limited_shadow: bool = False,
    min_az_sunrise: float = 120.0,
    max_az_sunrise: float = 290.0,
    time_to_sunrise: float = 3.0,
    sun_alt_limit: float | None = None,
    apply_cloud_mask: bool = False,
    cloud_limit: float = 1.5,
) -> list[bf.BaseBasisFunction]:
    """Basic standard mask basis functions.

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
        Should be set to the expected time needed to execute the
        observing block associated with the survey.
    apply_time_limited_shadow : `bool`, optional
        Flag for whether to apply the morning (time_to_sunrise) azimuth mask.
    min_az_sunrise : `float`, optional
        Minimum azimuth angle (in degrees) to observe during time period
        at the end of the night (during time_to_sunrise).
    max_az_sunrise: `float`, optional
        Maximum azimuth angle (in degrees) to observe during the time period
        at the end of the night (during time_to_sunrise).
    time_to_sunrise : `float`, optional
        Hours before daybreak (sun @ alt=0) to start the azimuth avoidance
        mask.
    sun_alt_limit : `float` or None, optional
        Maximum sun altitude (deg) required before proposing targets.
    apply_cloud_mask : `bool`, optional
        Flag for weather to add the cloud map mask, avoiding areas on sky
        with DREAM reported cloud values > `cloud_limit`.
    cloud_limit : `float`, optional
        The cloud extinction limit (in mag, from DREAM) required to trigger
        avoiding that area of the sky.
        In general, this should be enabled for surveys which only propose
        a few visits at a time, but disabled for surveys which
        generate long queues
        It will also be enabled at the QueueManager level.

    Returns
    -------
    mask_basis_functions : `list` [`BaseBasisFunction`]
        Mask basis functions should always be used with a weight of 0.
        The masked (np.nan or -np.inf) regions will remain masked,
        but the basis function values won't influence the reward.


    Notes
    -----
    Avoids the moon, bright planets, high wind, and
    areas on the sky out of bounds, using
    the MoonAvoidanceBasisFunction, PlanetMaskBasisFunction,
    MaskDirectWindBasisFunction, and the AltAzShadowMaskBasisFunction.

    If specified, sun_alt_limit will enforce only operating while the sun is
    below the sun_alt_limit (in degrees).

    If specified, a cloud avoidance mask will be added, using DREAM data.
    This is useful to add for surveys providing a short queue of observations,
    but is less useful for surveys with a long queue of observations, where it
    is better to implement (only) at the queue manager. Short vs. long is
    probably something like five minutes' worth, but will depend on how
    fast the clouds are moving.

    Optionally, adds the default AltAzShadowMaskTimeLimited basis function
    to avoid pointing toward sunrise during the last 3 hours of the night.
    This is now primarily used when the observatory is on generator power.
    """
    mask_bfs = []
    # Avoid the moon - too close to the moon will trip the REBs
    mask_bfs.append(
        bf.MoonAvoidanceBasisFunction(nside=nside, moon_distance=moon_distance)
    )
    # Avoid fast moving bright planets
    mask_bfs.append(bf.PlanetMaskBasisFunction(nside=nside))
    # Avoid the wind (mask only)
    mask_bfs.append(
        bf.MaskDirectWindBasisFunction(
            nside=nside, wind_speed_maximum=wind_speed_maximum
        )
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

    # If there are only short sequences, adding a cloud map here is good.
    if apply_cloud_mask:
        mask_bfs.append(
            bf.MaskCloudMapBasisFunction(nside=nside, extinction_limit=cloud_limit)
        )

    if apply_time_limited_shadow:
        # Only look away from the azimuth of the sun in the next day
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

    if sun_alt_limit is not None:
        mask_bfs.append(bf.SunAltLimitBasisFunction(alt_limit=sun_alt_limit))

    return mask_bfs


def standard_bf(
    nside: int = DEFAULT_NSIDE,
    bandname: str = "g",
    bandname2: str | None = "i",
    m5_weight: float = 6.0,
    fiducial_fwhm: float = 1.3,
    apply_cloud_extinction: bool = True,
    footprint_weight: float = 1.5,
    slewtime_weight: float = 3.0,
    stayband_weight: float = 3.0,
    footprints: Footprints | None = None,
    strict: bool = True,
    seeing_fwhm_max: float | npt.NDArray | None = None,
) -> list[tuple[bf.BaseBasisFunction, float]]:
    """Generate the standard basis functions that are shared by blob surveys

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
    fiducial_fwhm : `float`
        The fiducial FWHM for the M5Diff Basis function.
    apply_cloud_extinction : `bool`
        Turn the cloud extinction flag on (True) or off (False) in the M5Diff
        basis function.
    footprint_weight : `float`
        The weight on the survey footprint basis function.
    slewtime_weight : `float`
        The weight on the slewtime basis function.
    stayband_weight : `float`
        The weight on basis function that tries to stay avoid band changes.
    footprints : `rubin_scheduler.scheduler.utils.Footprints` object
        The desired footprints object. Default of None will work, but is likely
        not desirable.
    strict : `bool`
        If False, use BandChangeBasisFunction which rewards visits in the
        same bandpass as currently in-use.
        If True, use a StrictBandBasisFunction which rewards visits in the
        same bandpass as currently in-use, but also rewards visits in
        different filters if the moon rose/set or twilight ended/started,
        or if there was a large gap in observing.
    seeing_fwhm_max : `float` or `np.ndarray`
        Seeing limit to pass to the FootprintBasisFunction - visits
        with delivered image quality > seeing_fwhm_max will not be counted.

    Returns
    -------
    basis_functions_weights : `list`
        list of tuple pairs (basis function, weight) that is
        (rubin_scheduler.scheduler.BasisFunction object, float)

    """

    bfs = []

    if bandname2 is not None:
        bfs.append(
            (
                bf.M5DiffBasisFunction(
                    bandname=bandname,
                    nside=nside,
                    fiducial_FWHMEff=fiducial_fwhm,
                    apply_cloud_extinction=apply_cloud_extinction,
                ),
                m5_weight / 2.0,
            )
        )
        bfs.append(
            (
                bf.M5DiffBasisFunction(
                    bandname=bandname2,
                    nside=nside,
                    fiducial_FWHMEff=fiducial_fwhm,
                    apply_cloud_extinction=apply_cloud_extinction,
                ),
                m5_weight / 2.0,
            )
        )

    else:
        bfs.append(
            (
                bf.M5DiffBasisFunction(
                    bandname=bandname,
                    nside=nside,
                    fiducial_FWHMEff=fiducial_fwhm,
                    apply_cloud_extinction=apply_cloud_extinction,
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
                    seeing_fwhm_max=seeing_fwhm_max,
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
                    seeing_fwhm_max=seeing_fwhm_max,
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
                    seeing_fwhm_max=seeing_fwhm_max,
                ),
                footprint_weight,
            )
        )

    # Add a slewtime basis function. This probably should be bandpass = None.
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

    bandnames = [fn for fn in [bandname, bandname2] if fn is not None]
    bfs.append((bf.BandLoadedBasisFunction(bandnames=bandnames), 0))

    return bfs


def gen_template_surveys(
    footprints: Footprints,
    nside: int = DEFAULT_NSIDE,
    seeing_fwhm_max_zenith: dict = SEEING_FWHM_MAX_ZENITH_DEFAULT,
    median_cloud_limit: float = 1.5,
    band1s: list[str] = ["u", "g", "r", "i", "z", "y"],
    band2s: list[str] = ["u", "g", "r", "i", "z", "y"],
    dark_only: list[str] = ["u", "g"],
    ignore_obs: str | list[str] = ["DD", "twilight_near_sun"],
    camera_rot_limits: tuple[float, float] = CAMERA_ROT_LIMITS,
    exptime: float = EXPTIME,
    u_exptime: float = U_EXPTIME,
    n_obs_template: dict | None = None,
    pair_time: float = 25.0,
    area_required: float = 50.0,
    HA_min: float = 2.5,
    HA_max: float = 24 - 2.5,
    additional_area_limits: tuple[float] = (10.0,),
    extra_HA_mins: tuple[float] = (1.25,),
    extra_HA_maxes: tuple[float] = (24.0 - 1.25,),
    night_max: int = 365,
    m5_weight: float = 6.0,
    apply_cloud_extinction: bool = True,
    footprint_weight: float = 1.5,
    slewtime_weight: float = 3.0,
    stayband_weight: float = 3.0,
    science_program: str = SCIENCE_PROGRAM,
    blob_survey_params: dict | None = None,
    standard_mask_params: dict | None = None,
    pair_pad: float = 5.0,
) -> list[BlobSurvey]:
    """Surveys that are intended to acquire template visits in a convenient yet
    aggressive manner. Visits are aquired in pairs, with shorter than standard
    separation. There are additional constraints on delivered FWHM and
    overall cloud coverage (assuming that clouds in the sky imply
    that the sky background will be brighter due to scattering).

    Parameters
    ----------
    footprints : `rubin_scheduler.scheduler.utils.Footprints`
        The Footprints object for the Surveys.
    nside : `int`
        Nside for the surveys.
    seeing_fwhm_max_zenith : `dict` [ `str`, `float` ]
        The maximum (model delivered) seeing acceptable for the template
        surveys, at zenith. In arcseconds, as a dictionary one value per band..
        This limit is set to attempt to match the seeing limits
        used when building difference image coadds.
        When a template survey pair is composed of two different bandpasses,
        the minimum value between the pair will be used.
    median_cloud_limit : `float`
        The median cloud extinction limit (over the whole sky) before marking
        survey infeasible. Only applies if DREAM data available.
    band1s : `list` [`str`]
        The bandnames for the first band in a pair.
    band2s : `list` of `str`
        The band names for the second in the pair (None if unpaired).
    dark_only : `list` [ `str` ]
        The bands to only attempt during dark-time.
    ignore_obs : `str` or `list` [ `str` ]
        Strings to match within scheduler_note to flag observations to ignore.
    camera_rot_limits : `list` [ `float`, `float` ]
        Camera rotator limits (in degrees) for the dither rotation detailer.
    exptime : `float`
        The exposure time for grizy visits.
    u_exptime : `float`
        The exposure time for u band visits.
    n_obs_template : `dict` { `str` : `int` }
        Number of visits per bandpass before masking the Survey healpix.
        When this many visits are acquired, the pixel is considered "done".
    pair_time : `float`
        The time until the end of the first pass of the Blob.
        Since there is no second filter, this is the amount of time
        spent in the Blob.
    area_required : `float`
        The area required that needs templates, before the BlobSurvey will
        activate. Square degrees.
    HA_min : `float`
        The minimum HA to consider when considering template area.
    HA_max : `float`
        The maximum HA to consider when considering template area.
    additional_area_required : `list` [ `float` ]
        If the original area_required is not met, then allow a fallback
        to a second (smaller) area that will meet  the additional_HA
        specifications. This enables acquiring templates on small areas
        that did not fit into the original larger area coverage.
        Square degrees.
    extra_HA_mins : `float`
        When cleaning up the smaller area, additional_area_required,
        these small areas must fall within these HA requirements.
        Should be same length as additional_area_required.
    extra_HA_maxes : `float`
        When cleaning up the smaller area, additional_area_required,
        these small areas must fall within these HA requirements.
        Should be same length as additional_area_required.
    night_max : `int`
        The maximum number of nights after survey start to acquire templates.
    m5_weight : `float`
        The weight for the 5-sigma depth difference basis function.
    apply_cloud_extinction : `bool`
        Turn the cloud extinction flag on (True) or off (False) in the M5Diff
        basis function in the standard basis functions.
    footprint_weight : `float`
        The weight on the survey footprint basis function.
    slewtime_weight : `float`
        The weight on the slewtime basis function.
    stayband_weight : `float`
        The weight on basis function that tries to stay avoid band changes.
    science_program : `str`
        The science_program to use for visits from these surveys.
    blob_survey_params : `dict` or None
        A dictionary of additional kwargs to pass to the BlobSurvey.
        In particular, the times for typical slews, readtime, etc. are
        useful for setting the number of pointings to schedule within
        pair_time.
    standard_mask_params : `dict` or None
        A dictionary of additional kwargs to mass to the standard masks.
    pair_pad : `float`
        How much extra time (in minutes) to pad above the necessary pair time
        for shadow basis function.
    """

    if n_obs_template is None:
        n_obs_template = {"u": 5, "g": 5, "r": 5, "i": 5, "z": 5, "y": 5}

    if blob_survey_params is None:
        blob_survey_params = BLOB_SURVEY_PARAMS_DEFAULTS

    if standard_mask_params is None:
        standard_mask_params = STANDARD_MASK_DEFAULTS
        standard_mask_params["nside"] = nside
    else:
        standard_mask_params = copy.deepcopy(standard_mask_params)

    # Templates are acquired in pairs. Shadow mask accordingly.
    if band2s is None:
        shadow_minutes = pair_time + pair_pad
    else:
        shadow_minutes = pair_time * 2 + pair_pad
    if (
        "shadow_minutes" not in standard_mask_params
        or standard_mask_params["shadow_minutes"] < shadow_minutes
    ):
        standard_mask_params["shadow_minutes"] = shadow_minutes

    surveys = []

    for bandname, bandname2 in zip(band1s, band2s):

        # Set up Detailers for camera rotator and ordering in altitude.
        detailer_list = []
        detailer_list.append(
            detailers.CameraRotDetailer(
                min_rot=np.min(camera_rot_limits), max_rot=np.max(camera_rot_limits)
            )
        )
        detailer_list.append(detailers.CloseAltDetailer())
        # Make sure u band has different exposure time.
        detailer_list.append(
            detailers.BandNexp(bandname="u", nexp=1, exptime=u_exptime)
        )
        detailer_list.append(detailers.LabelRegionsAndDDFs())
        # Add extinction_limit detailer for cloud masking in queue_manager.
        detailer_list.append(
            detailers.ExtinctionLimitDetailer(
                extinction_limit=standard_mask_params["cloud_limit"]
            )
        )

        # For the bandpasses in use in this template survey,
        # find the seeing_fwhm_max using the minimum value for these bands.
        zenith_fwhm_max = min(
            seeing_fwhm_max_zenith[bandname], seeing_fwhm_max_zenith[bandname2]
        )
        # Calculate the seeing_fwhm_max array for the templates,
        # allowing for declination dependence.
        dec_fwhm_max = declination_dependent_fwhm(
            nside=nside, zenith_fwhm_limit=zenith_fwhm_max
        )

        # List to hold tuples of (basis_function_object, weight)
        bfs = []

        bfs.extend(
            standard_bf(
                nside,
                bandname=bandname,
                bandname2=bandname2,
                footprints=footprints,
                m5_weight=m5_weight,
                apply_cloud_extinction=apply_cloud_extinction,
                footprint_weight=footprint_weight,
                slewtime_weight=slewtime_weight,
                stayband_weight=stayband_weight,
                seeing_fwhm_max=dec_fwhm_max,
            )
        )
        # Add constraints around twilight.
        bfs.append((bf.TimeToTwilightBasisFunction(time_needed=pair_time), 0.0))
        bfs.append((bf.NotTwilightBasisFunction(), 0.0))
        # Push observations toward the meridian more strongly.
        bfs.append(
            (bf.RevHaMaskBasisFunction(ha_min=HA_min, ha_max=HA_max, nside=nside), 0.0)
        )

        # Need a once in night mask, so we spread epochs across many nights.
        bfs.append((bf.NInNightMaskBasisFunction(n_limit=1, nside=nside), 0.0))

        # If u or g, only when moon is down and outside of twilight.
        if bandname in dark_only:
            bfs.append((bf.NotTwilightBasisFunction(), 0.0))
            bfs.append((bf.MoonAltLimitBasisFunction(alt_limit=-5), 0.0))

        # Limit to first year
        bfs.append((bf.OnlyBeforeNightBasisFunction(night_max=night_max), 0.0))

        # Limit to only good seeing visits.
        bfs.append(
            (
                bf.MaskPoorSeeing(
                    bandname,
                    nside=nside,
                    seeing_fwhm_max=dec_fwhm_max,
                ),
                0,
            )
        )

        # Mask anything observed n_obs_template times.
        bfs.append(
            (
                bf.MaskAfterNObsSeeingBasisFunction(
                    nside=nside,
                    n_max=n_obs_template[bandname],
                    bandname=bandname,
                    seeing_fwhm_max=dec_fwhm_max,
                    reset_per_season=False,
                ),
                0.0,
            )
        )

        # Do not attempt if generally cloudy
        bfs.append(
            (
                bf.CloudedOutMapBasisFunction(median_cloud_limit=median_cloud_limit),
                0,
            )
        )

        # Add standard masks (including cloud mask)
        masks = standard_masks(**standard_mask_params)
        for m in masks:
            bfs.append((m, 0))

        # Unpack the basis functions and weights
        weights = [val[1] for val in bfs]
        basis_functions = [val[0] for val in bfs]

        # Additional masks - for the smaller additional_area_limits.
        additional_masks = []
        for ha_min, ha_max in zip(extra_HA_mins, extra_HA_maxes):
            additional_masks.append(
                bf.RevHaMaskBasisFunction(ha_min=ha_min, ha_max=ha_max, nside=nside)
            )

        survey_name = "templates %s%s" % (bandname, bandname2)
        observation_reason = f"template_blob_{bandname}{bandname2}_{pair_time:.1f}"

        surveys.append(
            BlobPairsSurvey(
                basis_functions,
                weights,
                bandname1=bandname,
                bandname2=bandname2,
                exptime=exptime,
                ideal_pair_time=pair_time,
                dither="call",
                survey_name=survey_name,
                science_program=science_program,
                observation_reason=observation_reason,
                ignore_obs=ignore_obs,
                nexp=1,
                detailers=detailer_list,
                area_required=area_required,
                additional_masks=additional_masks,
                additional_area_limits=additional_area_limits,
                note_block_size=True,
                **blob_survey_params,
            )
        )
    return surveys


def blob_for_long(
    footprints: Footprints,
    nside: int = DEFAULT_NSIDE,
    band1s: list[str] = ["g"],
    band2s: list[str] = ["i"],
    ignore_obs: str | list[str] = ["DD", "twilight_near_sun", "ToO"],
    camera_rot_limits: tuple[float, float] = CAMERA_ROT_LIMITS,
    exptime: float = EXPTIME,
    u_exptime: float = U_EXPTIME,
    pair_time: float = 33.0,
    HA_min: float = 12,
    HA_max: float = 24 - 3.5,
    m5_weight: float = 6.0,
    apply_cloud_extinction: bool = True,
    footprint_weight: float = 1.5,
    slewtime_weight: float = 3.0,
    stayband_weight: float = 3.0,
    night_pattern: tuple[bool, ...] = (
        True,
        True,
    ),
    time_after_twi: float = 30.0,
    blob_names: list[str] = [],
    scheduled_respect: float = 15.0,
    science_program: str = SCIENCE_PROGRAM,
    observation_reason: str | None = None,
    blob_survey_params: dict | None = None,
    standard_mask_params: dict | None = None,
    pair_pad: float = 5.0,
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
    u_exptime : `float`
        The exposure time for u band visits.
    pair_time : `float`
        The ideal time between pairs (minutes). Default 33.
    HA_min : `float`
        The minimum HA to consider when considering template area.
    HA_max : `float`
        The maximum HA to consider when considering template area.
    m5_weight : `float`
        The weight for the 5-sigma depth difference basis function.
    apply_cloud_extinction : `bool`
        Flag to include cloud extinction into M5Diff basis function.
    footprint_weight : `float`
        The weight on the survey footprint basis function.
    slewtime_weight : `float`
        The weight on the slewtime basis function.
    stayband_weight : `float`
        The weight on basis function that tries to stay avoid band changes.
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
    blob_survey_params : `dict` or None
        A dictionary of additional kwargs to pass to the BlobSurvey.
        In particular, the times for typical slews, readtime, etc. are
        useful for setting the number of pointings to schedule within
        pair_time.
    standard_mask_params : `dict` or None
        A dictionary of additional kwargs to mass to the standard masks.
    pair_pad : `float`
        How much extra time to pad above the necessary pair time
        for shadow basis function.
    """

    if blob_survey_params is None:
        blob_survey_params = BLOB_SURVEY_PARAMS_DEFAULTS
    if standard_mask_params is None:
        standard_mask_params = STANDARD_MASK_DEFAULTS
        standard_mask_params["nside"] = nside
    else:
        standard_mask_params = copy.deepcopy(standard_mask_params)

    # Calculate and apply shadow minutes to standard masks.
    if band2s is None:
        shadow_minutes = pair_time + pair_pad
    else:
        shadow_minutes = pair_time * 2 + pair_pad
    if (
        "shadow_minutes" not in standard_mask_params
        or standard_mask_params["shadow_minutes"] < shadow_minutes
    ):
        standard_mask_params["shadow_minutes"] = shadow_minutes

    surveys = []

    for bandname, bandname2 in zip(band1s, band2s):

        # Detailers.
        detailer_list = []
        detailer_list.append(
            detailers.CameraRotDetailer(
                min_rot=np.min(camera_rot_limits), max_rot=np.max(camera_rot_limits)
            )
        )
        detailer_list.append(detailers.CloseAltDetailer())
        detailer_list.append(
            detailers.BandNexp(bandname="u", nexp=1, exptime=u_exptime)
        )
        detailer_list.append(detailers.LabelRegionsAndDDFs())
        if standard_mask_params["apply_cloud_mask"]:
            # Add extinction_limit detailer for cloud masking in queue_manager.
            detailer_list.append(
                detailers.ExtinctionLimitDetailer(
                    extinction_limit=standard_mask_params["cloud_limit"]
                )
            )

        # List to hold tuples of (basis_function_object, weight)
        bfs = []

        bfs.extend(
            standard_bf(
                nside,
                bandname=bandname,
                bandname2=bandname2,
                m5_weight=m5_weight,
                apply_cloud_extinction=apply_cloud_extinction,
                footprint_weight=footprint_weight,
                slewtime_weight=slewtime_weight,
                stayband_weight=stayband_weight,
                footprints=footprints,
            )
        )

        # Make sure we respect scheduled observations
        bfs.append((bf.TimeToScheduledBasisFunction(time_needed=scheduled_respect), 0))
        # Add constraints around twilight
        bfs.append((bf.TimeToTwilightBasisFunction(time_needed=scheduled_respect), 0.0))
        bfs.append((bf.NotTwilightBasisFunction(), 0.0))
        bfs.append((bf.AfterEveningTwiBasisFunction(time_after=time_after_twi), 0.0))
        # Add HA constraints, so that first blob of triplet happens early
        bfs.append(
            (bf.HaMaskBasisFunction(ha_min=HA_min, ha_max=HA_max, nside=nside), 0.0)
        )
        # don't execute every night
        bfs.append((bf.NightModuloBasisFunction(night_pattern), 0.0))
        # only execute one blob per night
        bfs.append((bf.OnceInNightBasisFunction(notes=blob_names), 0))

        # Add standard masks
        masks = standard_masks(**standard_mask_params)
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
            observation_reason = f"triplet_pairs_{bandname}{bandname2}_{pair_time:.1f}"

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
                nexp=1,
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
    u_exptime: float = U_EXPTIME,
    pair_time: float = 33.0,
    night_pattern: tuple[bool, ...] = (True, False, False, False),
    gap_range: list[float] = [2, 7],
    HA_min: float = 12,
    HA_max: float = 24 - 3.5,
    time_after_twi: float = 120,
    m5_weight: float = 6.0,
    apply_cloud_extinction: bool = True,
    mask_cloud_limit: float = 1.5,
    footprint_weight: float = 1.5,
    slewtime_weight: float = 3.0,
    stayband_weight: float = 3.0,
    science_program: str = SCIENCE_PROGRAM,
    blob_survey_params: dict | None = None,
    standard_mask_params: dict | None = None,
    pair_pad: float = 5,
) -> list[LongGapSurvey]:
    """Generate long-gaps (triplets) surveys.

    Parameters
    -----------
    footprints : `rubin_scheduler.scheduler.utils.footprints.Footprints`
        The footprints to be used for the long-gaps surveys.
    nside : `int`
        The nside for the surveys.
    camera_rot_limits : `list` of `float`
        The limits to impose when rotationally dithering the camera (degrees).
    exptime : `float`
        The exposure time for grizy visits.
    u_exptime : `float`
        The exposure time for u band visits.
    pair_time : `float`
        The ideal time between pairs (minutes). Default 33.
    night_pattern : `list` [ `bool` ]
        Which nights to let the survey execute.
    gap_range : `list` [ `float` ]
        Range of times to attempt to gather pairs (hours).
    HA_min : `float`
        The hour angle limits passed to the initial blob scheduler. In hours.
    HA_max : `float`
        The hour angle limits passed to the initial blob scheduler.
    time_after_twi : `float`
        The time after evening twilight to attempt long gaps (minutes).
    m5_weight : `float`
        The weight for the 5-sigma depth difference basis function.
    apply_cloud_extinction : `bool`
        Flag to include cloud extinction into M5Diff basis function.
    mask_cloud_limit : `float`
        The extinction_limit to use for masking the survey, when
        cloud masking is active in the queue manager and survey.
    footprint_weight : `float`
        The weight on the survey footprint basis function.
    slewtime_weight : `float`
        The weight on the slewtime basis function.
    stayband_weight : `float`
        The weight on basis function that tries to stay avoid band changes.
    science_program : `str`
        The science_program to use for visits from these surveys.
    blob_survey_params : `dict` or None
        A dictionary of additional kwargs to pass to the BlobSurvey.
        In particular, the times for typical slews, readtime, etc. are
        useful for setting the number of pointings to schedule within
        pair_time.
    standard_mask_params : `dict` or None
        A dictionary of additional kwargs to mass to the standard masks.
    pair_pad : `float`
        How much extra time to pad above the necessary pair time
        for shadow basis function.
    """
    # Only copy the standard_mask_params here for the ScriptedSurvey.
    # The blob_for_long survey will copy/update on its own.
    if standard_mask_params is None:
        standard_mask_params_scripted: dict[str, Any] = {"nside": nside}
    else:
        standard_mask_params_scripted = copy.deepcopy(standard_mask_params)
    if (
        "shadow_minutes" not in standard_mask_params_scripted
        or standard_mask_params_scripted["shadow_minutes"] < pair_time
    ):
        standard_mask_params_scripted["shadow_minutes"] = pair_time

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
            u_exptime=u_exptime,
            pair_time=pair_time,
            nside=nside,
            band1s=[bandname1],
            band2s=[bandname2],
            night_pattern=night_pattern,
            time_after_twi=time_after_twi,
            HA_min=HA_min,
            HA_max=HA_max,
            m5_weight=m5_weight,
            footprint_weight=footprint_weight,
            slewtime_weight=slewtime_weight,
            stayband_weight=stayband_weight,
            blob_names=blob_names,
            science_program=science_program,
            blob_survey_params=blob_survey_params,
            standard_mask_params=standard_mask_params,
            pair_pad=pair_pad,
        )
        masks = standard_masks(**standard_mask_params_scripted)
        scripted = ScriptedSurvey(
            masks,
            nside=nside,
            ignore_obs=["blob", "DDF", "twi", "pair", "templates", "ToO"],
            science_program=science_program,
            detailers=[detailers.LabelRegionsAndDDFs()],
        )
        surveys.append(
            LongGapSurvey(blob[0], scripted, gap_range=gap_range, avoid_zenith=True)
        )

    return surveys


def gen_greedy_surveys(
    nside: int = DEFAULT_NSIDE,
    bands: list[str] = ["u", "g", "r", "i", "z", "y"],
    ignore_obs: list[str] = ["DD", "twilight_near_sun", "ToO"],
    camera_rot_limits: tuple[float, float] = CAMERA_ROT_LIMITS,
    exptime: float = EXPTIME,
    u_exptime: float = U_EXPTIME,
    shadow_minutes: float = 15.0,
    m5_weight: float = 3.0,
    apply_cloud_extinction: bool = True,
    footprint_weight: float = 0.75,
    slewtime_weight: float = 3.0,
    stayband_weight: float = 100.0,
    repeat_weight: float = -1.0,
    footprints: Footprints | None = None,
    science_program: str = SCIENCE_PROGRAM,
    standard_mask_params: dict | None = None,
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
    u_exptime : `float`
        The exposure time for u band visits.
    shadow_minutes : `float`
        Used to mask regions around zenith (minutes).
    m5_weight : `float`
        The weight for the 5-sigma depth difference basis function.
    apply_cloud_extinction : `bool`
        Flag to include cloud extinction into M5Diff basis function.
    footprint_weight : `float`
        The weight on the survey footprint basis function.
    slewtime_weight : `float`
        The weight on the slewtime basis function.
    stayband_weight : `float`
        The weight on basis function that tries to stay avoid band changes.
    repeat_weight : `float`
        Weight that enhances (+ values) or decreases (- values) the likelihood
        of revisiting the same pointing within a two-hour time gap.
    footprints : `rubin_scheduler.scheduler.utils.footprints.Footprints`
        The footprints to be used for the long-gaps surveys.
    science_program : `str`
        The science_program to use for visits from these surveys.
    standard_mask_params : `dict` or None
        A dictionary of additional kwargs to mass to the standard masks.
    """
    # Define the extra parameters that are used in the greedy survey. I
    # think these are fairly set, so no need to promote to utility func kwargs
    greed_survey_params = {
        "block_size": 1,
        "smoothing_kernel": None,
        "seed": 42,
        "camera": "LSST",
        "dither": "night",
    }
    if standard_mask_params is None:
        standard_mask_params = {"nside": nside}
    else:
        standard_mask_params = copy.deepcopy(standard_mask_params)
    if (
        "shadow_minutes" not in standard_mask_params
        or standard_mask_params["shadow_minutes"] < shadow_minutes
    ):
        standard_mask_params["shadow_minutes"] = shadow_minutes
    # The greedy survey could potentially add a cloud mask here.
    # However we should probably not, to avoid losing this backup tier.

    surveys = []
    detailer_list = [
        detailers.CameraRotDetailer(
            min_rot=np.min(camera_rot_limits), max_rot=np.max(camera_rot_limits)
        )
    ]
    detailer_list.append(detailers.LabelRegionsAndDDFs())
    # This is probably False (to allow greedy survey to always run).
    if standard_mask_params["apply_cloud_mask"]:
        # Add extinction_limit detailer for cloud masking in queue_manager.
        detailer_list.append(
            detailers.ExtinctionLimitDetailer(
                extinction_limit=standard_mask_params["cloud_limit"]
            )
        )

    if "u" in bands:
        detailer_list.append(
            detailers.BandNexp(bandname="u", nexp=1, exptime=u_exptime)
        )

    for bandname in bands:
        bfs = []
        bfs.extend(
            standard_bf(
                nside,
                bandname=bandname,
                bandname2=None,
                m5_weight=m5_weight,
                apply_cloud_extinction=apply_cloud_extinction,
                footprint_weight=footprint_weight,
                slewtime_weight=slewtime_weight,
                stayband_weight=stayband_weight,
                footprints=footprints,
                strict=False,
            )
        )

        # XXX-magic numbers
        bfs.append(
            (
                bf.VisitRepeatBasisFunction(
                    gap_min=0, gap_max=2 * 60.0, bandname=None, nside=nside, npairs=20
                ),
                repeat_weight,
            )
        )
        masks = standard_masks(**standard_mask_params)
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
                nexp=1,
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
    ignore_obs: str | list[str] = ["DD", "twilight_near_sun", "ToO"],
    camera_rot_limits: tuple[float, float] = CAMERA_ROT_LIMITS,
    exptime: float = EXPTIME,
    u_exptime: float = U_EXPTIME,
    pair_time: float = 33.0,
    max_pair_time: float = 40.0,
    m5_weight: float = 6.0,
    apply_cloud_extinction: bool = True,
    footprint_weight: float = 1.5,
    slewtime_weight: float = 3.0,
    stayband_weight: float = 3.0,
    repeat_weight: float = -20,
    good_seeing: dict = {"g": 3, "r": 3, "i": 3},
    good_seeing_weight: float = 3.0,
    seeing_fwhm_best: float = 0.8,
    m5_penalty_max: float = 0.5,
    survey_start: float = SURVEY_START_MJD,
    scheduled_respect: float = 15.0,
    science_program: str = SCIENCE_PROGRAM,
    blob_survey_params: dict | None = None,
    standard_mask_params: dict | None = None,
    pair_pad: float = 5.0,
) -> list[BlobSurvey]:
    """Generate surveys that take observations in blobs.

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
    u_exptime : `float`
        The exposure time for u band visits.
    pair_time : `float`
        The ideal time between pairs (minutes).
    max_pair_time : `float`
        The maximum time between pairs (minutes).
        The dynamic scaling won't scale blobs beyond this.
    m5_weight : `float`
        The weight for the 5-sigma depth difference basis function.
    apply_cloud_extinction : `bool`
        Flag to include cloud extinction into M5Diff basis function.
    footprint_weight : `float`
        The weight on the survey footprint basis function.
    slewtime_weight : `float`
        The weight on the slewtime basis function.
    stayband_weight : `float`
        The weight on basis function that tries to stay avoid band changes.
    repeat_weight : `float`
        Weight that enhances (+ values) or decreases (- values) the likelihood
        of revisiting the same pointing within a two-hour time gap.
    good_seeing : `dict`
        Number of good-seeing (<6") images per band for basis function
        directing good-atmospheric seeing into particular bands.
    good_seeing_weight : `float`
        The weight to place on getting good-seeing images in the
        good_seeing bands.
    seeing_fwhm_best : `float`
        The FWHM effective to use when counting visits for the "GoodSeeing"
        basis function (not for templates; for deblending and shapes).
    m5_penalty_max : `float`
        The maximum penalty in 5-sigma limiting depth to consider
        still good for the 'good seeing' images. (in mag).
    survey_start : `float`
        The mjd that the survey started (used for determining season for
        counting good seeing images within a season).
    scheduled_respect : `float`
        Ensure that blobs don't start within this many minutes of scheduled
        observations (from a ScriptedSurvey). Also used for start of twilight.
    science_program : `str`
        The science_program to use for visits from these surveys.
    blob_survey_params : `dict` or None
        A dictionary of additional kwargs to pass to the BlobSurvey.
        In particular, the times for typical slews, readtime, etc. are
        useful for setting the number of pointings to schedule within
        pair_time.
    standard_mask_params : `dict` or None
        A dictionary of additional kwargs to mass to the standard masks.
    pair_pad : `float`
        Padding to add to shadow mask. Default 5 (minutes).
    """
    if blob_survey_params is None:
        blob_survey_params = BLOB_SURVEY_PARAMS_DEFAULTS

    if standard_mask_params is None:
        standard_mask_params = {"nside": nside}
    else:
        standard_mask_params = copy.deepcopy(standard_mask_params)

    shadow_minutes = pair_time * 2 + pair_pad
    if (
        "shadow_minutes" not in standard_mask_params
        or standard_mask_params["shadow_minutes"] < shadow_minutes
    ):
        standard_mask_params["shadow_minutes"] = shadow_minutes

    surveys = []

    for bandname, bandname2 in zip(band1s, band2s):

        # Detailers.
        detailer_list = []
        detailer_list.append(
            detailers.CameraRotDetailer(
                min_rot=np.min(camera_rot_limits), max_rot=np.max(camera_rot_limits)
            )
        )
        detailer_list.append(detailers.CloseAltDetailer())
        if (bandname == "u") | (bandname2 == "u"):
            detailer_list.append(
                detailers.BandNexp(bandname="u", nexp=1, exptime=u_exptime)
            )
        detailer_list.append(detailers.FlushForSchedDetailer())
        detailer_list.append(detailers.LabelRegionsAndDDFs())
        if standard_mask_params["apply_cloud_mask"]:
            # Add extinction_limit detailer for cloud masking in queue_manager.
            detailer_list.append(
                detailers.ExtinctionLimitDetailer(
                    extinction_limit=standard_mask_params["cloud_limit"]
                )
            )

        # List to hold tuples of (basis_function_object, weight)
        bfs = []

        bfs.extend(
            standard_bf(
                nside,
                bandname=bandname,
                bandname2=bandname2,
                m5_weight=m5_weight,
                apply_cloud_extinction=apply_cloud_extinction,
                footprint_weight=footprint_weight,
                slewtime_weight=slewtime_weight,
                stayband_weight=stayband_weight,
                footprints=footprints,
            )
        )

        # Suppress revisits within 3 hours of the first pair.
        # Without this, we tend to repeat fields too quickly
        # but repeats after 3 hours could be useful, so let it expire then.
        bfs.append(
            (
                bf.VisitRepeatBasisFunction(
                    gap_min=0, gap_max=3 * 60.0, bandname=None, nside=nside, npairs=20
                ),
                repeat_weight,
            )
        )

        # Get some high quality FWHM images each season, for deblending.
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
                            seeing_fwhm_max=seeing_fwhm_best,
                            m5_penalty_max=m5_penalty_max,
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
                            seeing_fwhm_max=seeing_fwhm_best,
                            m5_penalty_max=m5_penalty_max,
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
                            seeing_fwhm_max=seeing_fwhm_best,
                            m5_penalty_max=m5_penalty_max,
                        ),
                        good_seeing_weight,
                    )
                )

        # Make sure we respect scheduled observations
        bfs.append((bf.TimeToScheduledBasisFunction(time_needed=scheduled_respect), 0))
        bfs.append((bf.TimeToTwilightBasisFunction(time_needed=scheduled_respect), 0.0))
        bfs.append((bf.NotTwilightBasisFunction(), 0.0))

        # Add standard masks
        masks = standard_masks(**standard_mask_params)
        for m in masks:
            bfs.append((m, 0))

        # unpack the basis functions and weights
        weights = [val[1] for val in bfs]
        basis_functions = [val[0] for val in bfs]

        # Set survey name
        if bandname2 is None:
            survey_name = "pair_%i, %s" % (pair_time, bandname)
        else:
            survey_name = "pair_%i, %s%s" % (pair_time, bandname, bandname2)
        if bandname2 is not None:
            detailer_list.append(detailers.TakeAsPairsDetailer(bandname=bandname2))

        observation_reason = f"pairs_{bandname}"
        if bandname2 is not None:
            observation_reason += f"{bandname2}"
        observation_reason += f"_{pair_time:.1f}_{max_pair_time:.1f}"

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
                nexp=1,
                detailers=detailer_list,
                note_block_size=True,
                max_pair_time=max_pair_time,
                **blob_survey_params,
            )
        )

    return surveys


def generate_twilight_near_sun(
    nside: int = DEFAULT_NSIDE,
    night_pattern: list[bool] | None = None,
    exptime: float = 15,
    ideal_pair_time: float = 5.0,
    max_airmass: float = 2.0,
    camera_rot_limits: tuple[float, float] = CAMERA_ROT_LIMITS,
    time_needed: float = 10.0,
    footprint_mask: npt.NDArray | float = 1,
    footprint_weight: float = 0.1,
    slewtime_weight: float = 3.0,
    stayband_weight: float = 3.0,
    band_dist_weight: float = 0.3,
    min_area: float | None = None,
    bands: str = "riz",
    n_repeat: int = 4,
    sun_alt_limit: float = -14.8,
    time_to_12deg: float = 25.0,
    slew_estimate: float = 4.5,
    max_elong: float = 60.0,
    ignore_obs: list[str] = ["DD", "pair", "long", "blob", "greedy", "template", "ToO"],
    science_program: str = SCIENCE_PROGRAM,
    standard_mask_params: dict | None = None,
) -> list[BlobSurvey]:
    """Generate a survey for observing NEO objects in twilight.

    This aims to cover ecliptic area within the LSST footprint, near the sun.

    Parameters
    ----------
    nside : `int`
        The HEALpix nside to use. Default to DEFAULT_NSIDE.
    night_pattern : `list` of `bool`
        A list of bools that set when the survey will be
        active. e.g., [True, False] for every-other night,
        [True, False, False] for every third night.
        Default None.
    exptime : `float`
        Exposure time of visits. Default 15.
    ideal_pair_time : `float`
        Ideal time between repeat visits (minutes).
        Default 5
    max_airmass : `float`
        Maximum airmass to attempt (unitless). Default 2.
    camera_rot_limits : `list` of `float`
        The camera rotation limits to use (degrees).
        Default [-80, 80].
    time_needed : `float`
        How much time should be available
        (e.g., before twilight ends) (minutes).
        Default 10
    footprint_mask : `np.ndarray` or `float`
        Mask to apply to the constructed ecliptic target mask (None).
        Keeps ecliptic area within this region.
        Default None
    footprint_weight : `float`
        Weight for footprint basis function. Default 0.1 (unitless).
    slewtime_weight : `float`
        Weight for slewtime basis function. Default 3 (unitless)
    stayband_weight : `float`
        Weight for staying in the same band basis function.
        Default 3 (unitless)
    band_dist_weight : `float`
        Weight for the basis function that tries to keep the distribution
        of visits across filters even.
    min_area : `float`
        The area that needs to be available before the survey will return
        observations (sq degrees). Default None.
    bands : `str`
        The bands to use, default 'riz'
    n_repeat : `int`
        The number of times a blob should be repeated, default 4.
    sun_alt_limit : `float`
        Do not start unless sun is higher than this limit (degrees).
        Default -14.8.
    time_to_sunrise : `float`
        Do not execute if time to sunrise is greater than (minutes).
        Default 25.
    slew_estimate : `float`
        An estimate of how long it takes to slew between
        neighboring fields (seconds). Default 4.5
    max_elong : `float`
        Maximum solar elongation to mini-allow survey to reach.
    ignore_obs : `str` or `list` [ `str` ]
        Ignore observations by surveys that include the given substring(s).
    standard_mask_params : `dict` or None
        A dictionary of additional kwargs to mass to the standard masks.
    """
    if standard_mask_params is None:
        standard_mask_params = {"nside": nside}
    else:
        standard_mask_params = copy.deepcopy(standard_mask_params)
    if (
        "shadow_minutes" not in standard_mask_params
        or standard_mask_params["shadow_minutes"] < ideal_pair_time * 4
    ):
        # pair_time * 4 (for the quad)
        standard_mask_params["shadow_minutes"] = ideal_pair_time * 4

    survey_name = "twilight_near_sun"
    observation_reason = "twilight_near_sun"
    footprint = ecliptic_area(nside=nside, mask=footprint_mask)
    constant_fp = ConstantFootprint(nside=nside)
    for bandname in bands:
        constant_fp.set_footprint(bandname, footprint)

    surveys = []
    for bandname in bands:
        detailer_list = []
        detailer_list.append(
            detailers.CameraRotDetailer(
                min_rot=np.min(camera_rot_limits), max_rot=np.max(camera_rot_limits)
            )
        )
        detailer_list.append(detailers.CloseAltDetailer())
        # Should put in a detailer so things start at lowest altitude
        detailer_list.append(
            detailers.TwilightTripleDetailer(
                slew_estimate=slew_estimate, n_repeat=n_repeat
            )
        )
        detailer_list.append(detailers.RandomBandDetailer(bands=bands))
        detailer_list.append(detailers.LabelRegionsAndDDFs())
        if standard_mask_params["apply_cloud_mask"]:
            # Add extinction_limit detailer for cloud masking in queue_manager.
            detailer_list.append(
                detailers.ExtinctionLimitDetailer(
                    extinction_limit=standard_mask_params["cloud_limit"]
                )
            )

        bfs = []

        bfs.append(
            (
                bf.FootprintBasisFunction(
                    bandname=bandname,
                    footprint=constant_fp,
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
        bfs.append((bf.StrictBandBasisFunction(bandname=bandname), stayband_weight))
        bfs.append((bf.BandDistBasisFunction(bandname=bandname), band_dist_weight))
        # Need a toward the sun, reward high airmass, with an
        # airmass cutoff basis function.
        bfs.append(
            (
                bf.NearSunHighAirmassBasisFunction(
                    nside=nside, max_airmass=max_airmass
                ),
                0,
            )
        )

        bfs.append((bf.BandLoadedBasisFunction(bandnames=bandname), 0))
        bfs.append(
            (
                bf.SolarElongationMaskBasisFunction(
                    min_elong=0.0, max_elong=max_elong, nside=nside
                ),
                0,
            )
        )

        bfs.append((bf.NightModuloBasisFunction(pattern=night_pattern), 0))
        # Do not attempt unless the sun is getting high
        bfs.append(
            (
                bf.CloseToTwilightBasisFunction(
                    max_sun_alt_limit=sun_alt_limit, max_time_to_12deg=time_to_12deg
                ),
                0,
            )
        )

        # Add standard masks
        masks = standard_masks(**standard_mask_params)
        for m in masks:
            bfs.append((m, 0))

        # unpack the basis functions and weights
        weights = [val[1] for val in bfs]
        basis_functions = [val[0] for val in bfs]

        # Set huge ideal pair time and use the detailer to cut down
        # the list of observations to fit twilight?
        surveys.append(
            BlobSurvey(
                basis_functions,
                weights,
                bandname1=bandname,
                bandname2=None,
                ideal_pair_time=ideal_pair_time,
                nside=nside,
                exptime=exptime,
                survey_name=survey_name,
                ignore_obs=ignore_obs,
                dither="night",
                nexp=1,
                detailers=detailer_list,
                twilight_scale=False,
                check_scheduled=False,
                area_required=min_area,
                observation_reason=observation_reason,
                science_program=science_program,
            )
        )
    return surveys
