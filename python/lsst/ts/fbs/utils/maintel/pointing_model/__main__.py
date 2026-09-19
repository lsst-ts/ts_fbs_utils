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

import argparse

import numpy as np
from astropy import units
from astropy.coordinates import Angle

from ... import get_data_dir
from .pointing_model import AltAzLimits, make_pointing_grid_hadec


def main() -> None:
    parser = argparse.ArgumentParser("make_pointing_model_grid")

    parser.add_argument(
        "--min_altitude",
        type=float,
        help="Minimum altitude, in degrees.",
        default=45.0,
    )
    parser.add_argument(
        "--max_altitude",
        type=float,
        help="Maximum altitude, in degrees.",
        default=80.0,
    )
    parser.add_argument(
        "--min_azimuth",
        type=float,
        help="Minimum azimuth, in degrees.",
        default=-180.0,
    )
    parser.add_argument(
        "--max_azimuth",
        type=float,
        help="Maximum azimuth, in degrees.",
        default=180,
    )
    parser.add_argument(
        "--dec_min",
        type=float,
        help="Minimum declination, in degrees.",
        default=-70,
    )
    parser.add_argument(
        "--dec_max",
        type=float,
        help="Maximum declination, in degrees.",
        default=10,
    )
    parser.add_argument(
        "--hour_angle_n_grid",
        type=int,
        nargs="+",
        help="A sequence of numbers that specify the number of poitings in RA for each declination.",
        default=[3, 4, 6, 8, 5, 5, 5],
    )
    parser.add_argument(
        "--block",
        type=str,
        help="Block name.",
        required=True,
    )

    args = parser.parse_args()

    altaz_limits = AltAzLimits(
        min_altitude=Angle(args.min_altitude, units.deg),
        max_altitude=Angle(args.max_altitude, units.deg),
        min_azimuth=Angle(args.min_azimuth, units.deg),
        max_azimuth=Angle(args.max_azimuth, units.deg),
    )

    altaz_grid = make_pointing_grid_hadec(
        altaz_limit=altaz_limits,
        dec_min=Angle(args.dec_min, units.deg),
        dec_max=Angle(args.dec_max, units.deg),
        hour_angle_n_grid=args.hour_angle_n_grid,
    )

    sort_data = np.argsort([az.value for (_, az) in altaz_grid])
    filename = get_data_dir() / f"pointing-tiles-{args.block.lower()}.txt"

    print(f"Writing pointing grid data to {filename}.")

    with open(filename, "w") as fp:
        fp.write("Name Alt Az\n")
        for i, index in enumerate(sort_data):
            (alt, az) = altaz_grid[index]
            fp.write(f"{args.block}_{i+1:03d} {alt.deg:6.2f} {az.deg:7.2f}\n")


if __name__ == "__main__":
    main()
