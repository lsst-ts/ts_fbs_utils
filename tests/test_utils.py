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

import os

from lsst.ts.fbs.utils import (
    get_data_dir,
    get_auxtel_tiles,
)


def test_get_data_dir() -> None:

    assert os.path.exists(get_data_dir())


def test_get_auxtel_tiles() -> None:

    auxtel_tiles = get_auxtel_tiles()
    assert "Survey" in auxtel_tiles.colnames
    assert "Name" in auxtel_tiles.colnames
    assert "RA" in auxtel_tiles.colnames
    assert "Dec" in auxtel_tiles.colnames
