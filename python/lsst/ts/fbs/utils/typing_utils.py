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

__all__ = ["AssertSurvey"]

import typing

from .target import Target
from .tiles import Tiles


class AssertSurvey(typing.Protocol):
    """Utility class for type checking.

    This class defines the interface for methods that asserts inputs for
    surveys.
    """

    def __call__(
        self,
        spec_targets: typing.List[Target],
        image_tiles: typing.List[Tiles],
        image_targest: typing.List[Target],
    ) -> None:
        pass
