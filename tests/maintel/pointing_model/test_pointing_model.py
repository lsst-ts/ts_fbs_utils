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

import astropy.units as u
import pytest
from astropy.coordinates import Angle, EarthLocation
from lsst.ts.fbs.utils.maintel.pointing_model import AltAzLimits, find_start_ha_for_dec


@pytest.mark.parametrize("dec", [Angle(-45, u.deg), Angle(-80, u.deg)])
def test_find_start_ha_for_dec(dec: Angle) -> None:
    location = EarthLocation.from_geodetic(
        lon=-70.747698 * u.deg, lat=-30.244728 * u.deg, height=2663.0 * u.m
    )
    altaz_limit = AltAzLimits(
        min_altitude=Angle(20, u.deg),
        max_altitude=Angle(85, u.deg),
        min_azimuth=Angle(-220, u.deg),
        max_azimuth=Angle(220, u.deg),
    )

    ha = find_start_ha_for_dec(altaz_limit=altaz_limit, dec=dec, location=location)
    print(f"{ha}")
    assert Angle(0, u.hourangle) > ha > Angle(-12, u.hourangle)
