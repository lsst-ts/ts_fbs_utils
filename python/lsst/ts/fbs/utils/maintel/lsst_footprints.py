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
from rubin_scheduler.scheduler.utils import (
    CurrentAreaMap,
    Footprint,
    Footprints,
    get_template_coverage,
    make_rolling_footprints,
)
from rubin_scheduler.site_models import Almanac
from rubin_scheduler.utils import DEFAULT_NSIDE, SURVEY_START_MJD


def get_footprints(
    nside: int = DEFAULT_NSIDE,
    survey_start_mjd: float = SURVEY_START_MJD,
    bandpasses: tuple[str, ...] = ("u", "g", "r", "i", "z", "y"),
    roll_n_slice: int = 2,
    roll_scale: float = 0.9,
    roll_n_cycles: int = 3,
    roll_uniform: bool = True,
) -> tuple[Footprints, Footprint, npt.NDArray]:
    """Generate footprints for LSST surveys.

    Parameters
    -----------
    nside
        The NSIDE for the resolution of the footprint.
    survey_start_mjd
        The MJD of the survey start. This influences the footprint
        returned at a given time, as there is a slope in the footprint
        value over the observing season.
    bandpasses
        The list of bandpasses to include in the survey footprint.
        Deactivating or restricting bandpasses should be done here.
    roll_n_slice
        The number of slices to split the footprint into for the rolling
        cadence. 2 slices actually means 4 regions on the sky.
    roll_scale
        The strength to pass for the rolling cadence.
    roll_n_cycles
        The number of rolling cycles to include in the ten year survey
        footprint. One cycle = high activity + low activity cycle.
    roll_uniform
        Whether to include uniform rolling seasons in the rolling footprint
        generation.

    Returns
    -------
    lsst_footprints, template_footprint, footprints_mask :
        `Footprints`, `Footprint`, `np.NDArray`
        The standard LSST footprints to use for most surveys, that
        includes the WFD/NES/SCP/GP etc plus rolling cadence;
        the template footprint that already has the known templated
        area cut out;
        and a mask to constrain the ToO and near-sun twilight microsurvey.
    """

    # Generate footprint over the sky
    sky = CurrentAreaMap(nside=nside)
    footprints_hp_array, labels = sky.return_maps()
    # Repackage to dictionary because some downstream things like dictionaries.
    footprints_hp = {}
    for key in footprints_hp_array.dtype.names:
        # Don't mask based on bandpasses yet.
        footprints_hp[key] = footprints_hp_array[key]

    # Set up a mask to contain some surveys within this region.
    footprint_mask = footprints_hp["r"] * 0
    footprint_mask[np.where(footprints_hp["r"] > 0)] = 1

    # Identify pixels for rolling
    roll_indx = np.where((labels == "lowdust") | (labels == "virgo"))[0]
    roll_footprint = footprints_hp_array["r"] * 0
    roll_footprint[roll_indx] = 1

    # Now that we know where the rolling footprint is, we can mask if needed.
    for key in footprints_hp_array.dtype.names:
        if key not in bandpasses:
            # Zero out any footprints that are not in the bandpasses list.
            footprints_hp[key] = footprints_hp_array[key] * 0.0

    # Use the Almanac to find the position of the sun at the start of survey.
    almanac = Almanac(mjd_start=survey_start_mjd)
    sun_moon_info = almanac.get_sun_moon_positions(survey_start_mjd)
    sun_ra_start = sun_moon_info["sun_RA"].copy()

    # Define the rolling footprint
    footprints = make_rolling_footprints(
        fp_hp=footprints_hp,
        mjd_start=survey_start_mjd,
        sun_ra_start=sun_ra_start,
        nslice=roll_n_slice,
        scale=roll_scale,
        nside=nside,
        wfd_indx=roll_indx,
        order_roll=1,
        n_cycles=roll_n_cycles,
        uniform=roll_uniform,
    )

    # Create template footprint.
    # Similar to rolling footprint but tracks visits separately
    # (only good seeing visits) and no rolling.
    template_fp = Footprint(survey_start_mjd, sun_ra_start, nside=nside)
    # Read already-acquired templates from disk and remove from template fp.
    known_templates = get_template_coverage(nside=nside)
    # Combine goal footprint and known footprint
    for key in footprints_hp_array.dtype.names:
        if key not in bandpasses:
            tmp_fp = footprints_hp_array[key] * 0.0
        else:
            # Make template footprint flat across the whole known footprint.
            tmp_fp = np.where(footprints_hp_array[key] > 0, 1.0, 0.0)
            # Remove known templated areas
            tmp_fp = np.where(known_templates[key] > 0, 0.0, tmp_fp)
        template_fp.set_footprint(key, tmp_fp)

    # Set up a mask to contain ToO and neomicro surveys within LSST footprint.
    r_indx = footprints.bands["r"]
    footprint_mask = np.where(footprints.footprints[r_indx] > 0, 1.0, 0.0)

    return footprints, template_fp, footprint_mask
