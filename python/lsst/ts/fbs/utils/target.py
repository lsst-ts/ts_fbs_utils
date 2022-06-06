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

__all__ = ["Target"]

import typing

from astropy.coordinates import Angle
from dataclasses import dataclass


@dataclass
class Target:
    """Define a single target for the feature based scheduler.

    Attributes
    ----------
    target_name : `str`
        Name of the target.
    survey_name : `str`
        Name of the survey.
    ra : `Angle`
        Right ascention of the target as an astropy angle. Units are inferred.
    dec : `Angle`
        Declination of the target as an astropy angle. Units are inferred.
    hour_angle_limit : `list` of `tuple` of (`float`, `float`)
        List of tuples witth hour angle limits (min, max) in hour angle.
    reward_value : `float`
        Reward value for the target.
    filters: `list` of `str`
        List of filters to observe.
    visit_gap : `float`
        How long to wait wait until observing this target again after one
        observation, in minutes.
    exptime : `float`
        Total exposure time, in seconds.
    nexp : `int`
        Number of exposures.

    Notes
    -----
    The exposure time for each exposure is `exptime`/`nexp`.
    """

    target_name: str
    survey_name: str
    ra: Angle
    dec: Angle
    hour_angle_limit: typing.List[typing.Tuple[float, float]]
    reward_value: float
    filters: typing.List[str]
    visit_gap: float
    exptime: float
    nexp: int
