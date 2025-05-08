# This file is part of ts_fbs_utils
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
    "get_data_dir",
    "get_auxtel_tiles",
    "get_maintel_tiles",
    "get_pointing_model_grid_data",
]

import os
import pathlib

from astropy.io import ascii
from astropy.table import Table


def get_data_dir() -> pathlib.Path:
    """Return the data directory of this package."""
    if "LSST_TS_FBS_UTILS_DATA_PATH" in os.environ:
        data_path = pathlib.Path(os.environ["LSST_TS_FBS_UTILS_DATA_PATH"])
        if data_path.exists():
            return data_path.resolve().parents[0]

    return pathlib.Path(__file__).resolve().parents[0] / "data"


def get_auxtel_tiles() -> Table:
    """Get auxtel tiles.

    Returns
    -------
    `Table`
        AuxTel tiles as astropy table.
    """
    return ascii.read(
        str(get_data_dir() / "auxtel_tiles.txt"),
        format="basic",
        fast_reader=False,
        guess=False,
    )


def get_maintel_tiles() -> Table:
    """Get maintel tiles.

    Returns
    -------
    `Table`
        MainTel tiles as astropy table.
    """
    return ascii.read(
        str(get_data_dir() / "maintel_tiles.txt"),
        format="basic",
        fast_reader=False,
        guess=False,
    )


def get_pointing_model_grid_data(science_program: str) -> Table:
    """Get pointing model tiles.

    Parameters
    ----------
    science_program : `str`
        Name of the science program.

    Returns
    -------
    `Table`
        Pointing tiles as astropy tables.
    """
    return ascii.read(
        str(get_data_dir() / f"pointing-tiles-{science_program.lower()}.txt"),
        format="basic",
        fast_reader=False,
        guess=False,
    )
