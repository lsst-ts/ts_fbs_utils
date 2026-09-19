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

from dataclasses import dataclass

import astropy.units as u
import numpy as np
from astropy.coordinates import ICRS, AltAz, Angle, EarthLocation
from lsst.ts.utils import astropy_time_from_tai_unix, current_tai


@dataclass
class AltAzLimits:
    """
    A data class to store Altitude-Azimuth (AltAz)
    limits using `astropy.coordinates.Angle`.

    Parameters
    ----------
    min_altitude : Angle
        Minimum altitude. Must be between -90 and +90 degrees.
    max_altitude : Angle
        Maximum altitude. Must be between -90 and +90 degrees.
    min_azimuth : Angle
        Minimum azimuth. Typically between 0 and 360 degrees.
    max_azimuth : Angle
        Maximum azimuth. Typically between 0 and 360 degrees.

    Attributes
    ----------
    min_altitude : Angle
        Minimum altitude angle.
    max_altitude : Angle
        Maximum altitude angle.
    min_azimuth : Angle
        Minimum azimuth angle.
    max_azimuth : Angle
        Maximum azimuth angle.

    Raises
    ------
    ValueError
        If altitude values are out of range (-90 to +90 degrees),
        if `min_altitude` is greater than `max_altitude` or if
        `min_azimuth` is greater than `max_azimuth`.
    """

    min_altitude: Angle
    max_altitude: Angle
    min_azimuth: Angle
    max_azimuth: Angle

    def __post_init__(self) -> None:
        """Validates the angle limits after initialization."""
        self.min_altitude = self.min_altitude.to(u.deg)
        self.max_altitude = self.max_altitude.to(u.deg)
        self.min_azimuth = self.min_azimuth.to(u.deg)
        self.max_azimuth = self.max_azimuth.to(u.deg)

        if not (-90 * u.deg <= self.min_altitude <= 90 * u.deg):
            raise ValueError("min_altitude must be between -90 and +90 degrees.")
        if not (-90 * u.deg <= self.max_altitude <= 90 * u.deg):
            raise ValueError("max_altitude must be between -90 and +90 degrees.")
        if self.min_altitude > self.max_altitude:
            raise ValueError("min_altitude cannot be greater than max_altitude.")

        if self.min_azimuth > self.max_azimuth:
            raise ValueError("min_azimuth cannot be greater than max_azimuth.")

    def to_radians(self) -> None:
        """Convert altitude and azimuth limits to radians.

        Notes
        -----
        This method modifies the `min_altitude`, `max_altitude`, `min_azimuth`,
        and `max_azimuth` attributes in place.
        """
        self.min_altitude = self.min_altitude.to(u.rad)
        self.max_altitude = self.max_altitude.to(u.rad)
        self.min_azimuth = self.min_azimuth.to(u.rad)
        self.max_azimuth = self.max_azimuth.to(u.rad)

    def to_degrees(self) -> None:
        """Convert altitude and azimuth limits to degrees.

        Notes
        -----
        This method modifies the `min_altitude`, `max_altitude`, `min_azimuth`,
        and `max_azimuth` attributes in place.
        """
        self.min_altitude = self.min_altitude.to(u.deg)
        self.max_altitude = self.max_altitude.to(u.deg)
        self.min_azimuth = self.min_azimuth.to(u.deg)
        self.max_azimuth = self.max_azimuth.to(u.deg)

    def in_range(self, position: AltAz) -> bool:
        """Check if the input position is inside the limits.

        Parameters
        ----------
        position : `AltAz`
            Position to check.

        Returns
        -------
        `bool`
            If position is inside the limit.
        """
        return (
            self.min_altitude <= position.alt <= self.max_altitude
            and self.min_azimuth <= position.az <= self.max_azimuth
        )

    def in_range_alt_az(self, alt: Angle, az: Angle) -> bool:
        """Check if the input alt/az values are in range.

        Parameters
        ----------
        alt : `Angle`
            Altitude.
        az : `Angle`
            Azimuth.

        Returns
        -------
        `bool`
            If position is inside the limit.
        """
        return (
            self.min_altitude <= alt <= self.max_altitude
            and self.min_azimuth <= az <= self.max_azimuth
        )

    def in_alt_range(self, alt: Angle) -> bool:
        """Check if the input alt is in range.

        Parameters
        ----------
        alt : `Angle`
            Altitude.

        Returns
        -------
        `bool`
            If position is inside the limit.
        """
        return self.min_altitude <= alt <= self.max_altitude

    def approx_in_alt_range(self, alt: Angle, tol: Angle = Angle(0.1, u.deg)) -> bool:
        """Check if the input alt is in range.

        Parameters
        ----------
        alt : `Angle`
            Altitude.

        Returns
        -------
        `bool`
            If position is inside the limit.
        """
        return (
            self.min_altitude <= alt <= self.max_altitude
            or np.isclose(alt.deg, self.min_altitude.deg, atol=tol.deg)
            or np.isclose(alt.deg, self.max_altitude.deg, atol=tol.deg)
        )

    def __str__(self) -> str:
        """Return a string representation of the AltAzLimits instance.

        Returns
        -------
        str
            A string showing the altitude and azimuth limits
            in their current units.
        """
        return (
            f"AltAzLimits("
            f"min_altitude={self.min_altitude}, "
            f"max_altitude={self.max_altitude}, "
            f"min_azimuth={self.min_azimuth}, "
            f"max_azimuth={self.max_azimuth})"
        )


def make_pointing_grid_azel(
    altaz_limit: AltAzLimits,
    n_alt_min: int,
    n_alt_max: int,
    n_alt: int,
) -> list[AltAz]:
    """Build a pointing grid using azel as a reference frame.

    Parameters
    ----------
    altaz_limit : `AltAzLimits`
        Limits in alt/az for the grid.
    n_alt_min : `int`
        Minimum number of pointings in elevation.
    n_alt_max : `int`
        Maximum number of pointings in elevation.
    n_alt : `int`
        How many strips in elevation to construct?

    Returns
    -------
    `list`[`AltAz`]
        List of AltAz coordinates of the grid.
    """

    alt_grid = np.linspace(altaz_limit.min_altitude, altaz_limit.max_altitude, n_alt)
    altaz_grid = []

    for alt in alt_grid:
        n_az = int(
            np.ceil(
                np.max(
                    [
                        n_alt_max * np.cos(alt.rad - altaz_limit.min_altitude.rad),
                        n_alt_min,
                    ]
                )
            )
        )
        altaz_grid += [
            AltAz(alt=alt, az=az)
            for az in np.linspace(
                altaz_limit.min_azimuth, altaz_limit.max_altitude, n_az + 1
            )
            if altaz_limit.min_azimuth < az
        ]

    sort_index = np.argsort([altaz.az.value for altaz in altaz_grid])

    return [altaz_grid[index] for index in sort_index]


def make_pointing_grid_hadec(
    altaz_limit: AltAzLimits,
    dec_min: Angle,
    dec_max: Angle,
    hour_angle_n_grid: list[int],
) -> list[tuple[Angle, Angle]]:
    """Build a grid using hour angle/dec as a reference frame.

    Parameters
    ----------

    Returns
    -------
    `list`[`AltAz`]
        List of AltAz coordinates of the grid.
    """

    altaz_grid: list[tuple[Angle, Angle]] = []

    reverse_ha = False
    dec_values = np.linspace(
        dec_min,
        dec_max,
        len(hour_angle_n_grid),
    )

    location = EarthLocation.from_geodetic(
        lon=-70.747698 * u.deg, lat=-30.244728 * u.deg, height=2663.0 * u.m
    )

    for dec, n_ha in zip(dec_values, hour_angle_n_grid):
        ha_start = find_start_ha_for_dec(
            dec=dec,
            altaz_limit=altaz_limit,
            location=location,
        )
        ha_values = np.linspace(ha_start, -ha_start, n_ha)
        if reverse_ha:
            ha_values = ha_values[::-1]

        reverse_ha = not reverse_ha

        for ha in ha_values:
            altaz = azel_from_radec(location=location, ha=ha, dec=dec)
            alt = altaz.alt
            az = altaz.az.wrap_at(180.0 * u.deg)
            if altaz_limit.approx_in_alt_range(alt):
                altaz_grid.append((alt, az))

    return altaz_grid


def find_start_ha_for_dec(
    dec: Angle,
    altaz_limit: AltAzLimits,
    location: EarthLocation,
    threshold: Angle = Angle(0.1, u.deg),
    max_iter: int = 20,
) -> Angle:
    """Find the minimum hour angle that is inside the provided alt/az limits
    for a given declination.

    Parameters
    ----------
    dec : `Angle`
        Declination.
    altaz_limit : `AltAzLimits`
        The alt/az limits.
    threshold : `Angle`
        The threshold to find the start ha.
    max_iter : int
        Maximum number of iterations.

    Returns
    -------
    ha : `Angle`
        The minimum hour angle that is inside the alt/az limit.

    Raises
    ------
    ValueError
        If it cannot find an hour angle that satisfy the conditions.
    """

    max_hour_angle = Angle(0, u.hourangle)
    min_hour_angle = Angle(-12, u.hourangle)

    target_alt = altaz_limit.min_altitude

    max_altaz = azel_from_radec(location, max_hour_angle, dec)
    min_altaz = azel_from_radec(location, min_hour_angle, dec)

    if max_altaz.alt < target_alt:
        raise ValueError(
            f"Maximum altitude {max_altaz.alt.value:.2f} below target altitude {target_alt.value}."
        )

    if min_altaz.alt > target_alt:
        return min_hour_angle + threshold

    for i in range(max_iter):
        ha = (max_hour_angle + min_hour_angle) / 2.0
        current_altaz = azel_from_radec(location, ha, dec)
        if abs(current_altaz.alt - target_alt) < threshold:
            return ha
        elif current_altaz.alt > target_alt:
            max_hour_angle = ha
        else:
            min_hour_angle = ha
    else:
        raise ValueError("Could not find suitable hour angle.")


def azel_from_radec(location: EarthLocation, ha: Angle, dec: Angle) -> AltAz:
    """Calculate the AzEl position from an hour angle and dec."""

    time = astropy_time_from_tai_unix(current_tai())
    time.location = location
    sidereal_time = time.sidereal_time("mean")
    return ICRS(sidereal_time - ha, dec).transform_to(
        AltAz(location=location, obstime=time)
    )
