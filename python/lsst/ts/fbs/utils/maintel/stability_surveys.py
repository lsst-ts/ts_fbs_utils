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
    "gen_az_el_rot_stability_survey",
]

from typing import TypedDict

import numpy as np
from astropy.time import Time
from rubin_scheduler.scheduler.basis_functions import (
    BalanceVisits,
    RewardNObsSequence,
)
from rubin_scheduler.scheduler.detailers import (
    AltAz2RaDecDetailer,
    CopyValueDetailer,
    Rottep2RotspDesiredDetailer,
)
from rubin_scheduler.scheduler.surveys import FieldAltAzSurvey, FieldSurvey
from rubin_scheduler.utils import DEFAULT_NSIDE

from .lsst_surveys import EXPTIME, safety_masks


class StabilityTarget(TypedDict):
    note: str
    name: str
    alt: float
    az: float
    rotTelPos: float


def gen_az_el_rot_stability_survey(
    az_values: list[float],
    el_values: list[float],
    rot_values: list[float],
    science_program: str,
    observation_reason: str = "fbs driven aos stability test",
    sequence: list[str] = ["i"],
    nvis_per_cycle: int = 100,
    exptimes: dict | float = EXPTIME,
    nside: int = DEFAULT_NSIDE,
    safety_mask_params: dict | None = None,
) -> list[FieldSurvey]:
    """Generate a set of stability surveys, designed to execute a series of
    images at predefined Az El Rot positions. The result will be a list of
    FieldAltAzSurvey for each position combination provided.

    Parameters
    ----------
    az_values : list of `float`
        Azimuth values in degrees.
    el_values : list of `float`
        Elevation values in degrees.
    rot_values : list of `float`
        Rotator values (rotTelPos) in degrees.
    science_program : `str`
        Name of the science program for the survey.
    observation_reason : `str`
        Str to use for the observation reason for the survey.
    sequence : list[str], default ["i"]
        Band names to use in the sequence. Note that 1 image will be taken per
        band listed in sequences. It is recommended to only use 1 filter in the
        sequence to prevent the coordinate drifting too much between successive
        targets.
    nvis_per_cycle : `dict` or `int`
        Number of visits per position cycle. Total number of exposures at each
        position cycle is len(sequence) * nvis_per_cycle.
    exptimes : `dict` or `float`
        Exposure times for each band. Passed to FieldSurvey.
        If exptimes is a float, this is applied to each band in `sequence`.
        Default is 30s, this should be updated if a different exposure time is
        needed (e.g. 38s for u band).
    nside : `int`
        The HEALpix nside for the survey, used for basis functions.
    safety_mask_params : `dict`
        A dictionary of kwargs to pass to the standard safety masks.


    Returns
    -------
    surveys : `list` [`FieldSurvey`]
        A list of `FieldAltAzSurvey` objects, one per Az/El/Rot combination.
    """

    if safety_mask_params is None:
        safety_mask_params = {
            "nside": nside,
            "wind_speed_maximum": None,
            "shadow_minutes": 0,
            "apply_time_limited_shadow": False,
            "time_to_sunrise": 3.0,
            "min_az_sunrise": 144,
            "max_az_sunrise": 255,
        }

    day_obs = int(
        Time(int(Time.now().mjd - 0.5), format="mjd", scale="utc")
        .iso[0:10]
        .replace("-", "")
    )
    scheduler_root = f"{science_program} {day_obs}"

    target_dict: dict[str, StabilityTarget] = {}

    for alt in el_values:
        for az in az_values:
            for rotTelPos in rot_values:
                name = f"alt:{alt:.1f} az:{az:.1f} rotTel:{rotTelPos:.0f}"
                target_dict[name] = {
                    "note": f"{scheduler_root} {name}",
                    "name": f"{name}",
                    "alt": alt,
                    "az": az,
                    "rotTelPos": rotTelPos,
                }

    n_pointings = len(target_dict)

    detailers = [
        AltAz2RaDecDetailer(),
        Rottep2RotspDesiredDetailer(),
        CopyValueDetailer(source="rotSkyPos_desired", destination="rotSkyPos"),
    ]

    safety_masks_basis_functions = safety_masks(**safety_mask_params)

    survey_lists = []
    for target in target_dict:
        tt = target_dict[target]
        survey = FieldAltAzSurvey(
            basis_functions=[
                RewardNObsSequence(
                    n_obs_survey=nvis_per_cycle,
                    note_survey=tt["note"],
                ),
                BalanceVisits(
                    nobs_reference=nvis_per_cycle * n_pointings,
                    note_survey=tt["note"],
                    note_interest=f"{scheduler_root}",
                    nside=nside,
                ),
            ]
            + safety_masks_basis_functions,
            alt=tt["alt"],
            az=tt["az"],
            sequence=sequence,
            nvisits=1,
            exptimes=exptimes,
            ignore_obs=None,
            survey_name=tt["note"],
            target_name=tt["name"],
            science_program=science_program,
            observation_reason=observation_reason,
            scheduler_note=tt["note"],
            nside=nside,
            flush_pad=30.0,
            detailers=detailers,
        )
        survey.observations["rotTelPos"] = np.radians(tt["rotTelPos"])
        survey_lists.append(survey)

    return survey_lists
