################################################################################
# Copyright (C) 2014, 2015 GenAP, McGill University and Genome Quebec Innovation Centre
#
# This file is part of MUGQIC Pipelines.
#
# MUGQIC Pipelines is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MUGQIC Pipelines is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with MUGQIC Pipelines.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

#!/usr/bin/env python

# Python Standard Modules

# MUGQIC Modules
from core.job import Job

def mkdir(folder, remove=False):
    return Job(
        [],
        [folder],
        [],
        command="""\
mkdir -p {directory}""".format(
            directory=folder
        ),
        removable_files=[folder] if remove else []
    )

def ln(target_file, link, sleep=None):
    return Job(
        [target_file],
        [link],
        [],
        command="""\
ln -s -f \\
  {target_file} \\
  {link} \\
  {sleep}""".format(
            source=target_file,
            destination=link,
            sleep="&& sleep "+str(sleep) if sleep else ""
        ),
        removable_files=[link]
    )

def mv(source, target):
    return Job(
        [source],
        [target],
        [],
        command="""\
mv {source} {dest}""".format(
            source=source,
            dest=target
        )
    )

