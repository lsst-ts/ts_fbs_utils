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


import pathlib

import yaml


def get_data_dir() -> pathlib.Path:
    """Return the data directory of this package."""
    return pathlib.Path(__file__).resolve().parents[0]


def get_comcam_sv_targets() -> dict:
    infile = str(get_data_dir() / "field_survey_centers.yaml")
    with open(infile) as stream:
        try:
            targets_dict = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)

    return targets_dict["comcam_sv_targets"]
