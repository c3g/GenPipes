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

# Python Standard Modules
import os

# MUGQIC Modules
from core.job import Job

def mkdir(folder, remove=False, output_dependency=None):
    if output_dependency is not None:
        outputs=output_dependency
    else:
        outputs=[folder]
    return Job(
        [],
        outputs,
        [],
        command="""\
mkdir -p {directory} && \\
touch {directory}""".format(
            directory=folder
        ),
        removable_files=[folder] if remove else []
    )

def chgdir(folder):
    return Job(
        [],
        [folder],
        command="""\
cd {directory}""".format(
            directory=folder
        ),
    )

def ln(
    target_file,
    link,
    out_dir=None
    ):
    
    return Job(
        [target_file],
        [link],
        command="""\
ln -s -f \\
  {target_file} \\
  {link}""".format(
            target_file=os.path.join(out_dir, target_file) if out_dir else target_file,
            link=os.path.join(out_dir, link) if out_dir else link
        ),
        removable_files=[link]
    )

def mv(source, target, force=False, input_dependency=None, output_dependency=None):
    if input_dependency is not None:
        inputs=input_dependency
    else:
        inputs=[source]
    if output_dependency is not None:
        outputs=output_dependency
    else:
        outputs=[target]
    return Job(
        inputs,
        outputs,
        [],
        command="""\
mv {force}{source} \\
   {dest}""".format(
            force="-f " if force else "",
            source=source,
            dest=target
        )
    )

def cp(source, target, recursive=False, update=False, input_dependency=None, output_dependency=None):
    if input_dependency is not None:
        inputs=input_dependency
    else:
        inputs=[source]
    if output_dependency is not None:
        outputs=output_dependency
    else:
        outputs=[target]
    return Job(
        inputs,
        outputs,
        [],
        command="""\
cp {rec}{upd}'{source}' {dest}""".format(
            source=source,
            dest=target,
            rec="-r " if recursive else "",
            upd="-u " if update else ""
        )
    )

def rm(target, recursive=False, remove=False, force=False, input_dependency=None, output_dependency=[]):
    if input_dependency is not None:
        inputs=input_dependency
    else:
        inputs=[target]
    return Job(
        inputs,
        output_dependency,
        [],
        command="""\
rm {rec} {force} {target}""".format(
            target=target,
            rec="-r" if recursive else "",
            force="-f" if force else ""
        )
    )

def touch(target, input_dependency=None):
    if input_dependency is not None:
        inputs=input_dependency
    else:
        inputs=[target]
    return Job(
        inputs,
        [],
        [],
        command="""\
touch {target}""".format(
            target=target
        )
    )

def md5sum(
    inputs,
    output,
    check=False
    ):

    if not isinstance(inputs, list):
        inputs = [inputs]

    return Job(
        inputs,
        [output],
        command="""\
md5sum {check}{input}{output}""".format(
            check="-c " if check else "",
            input=" ".join(["'"+os.path.abspath(input)+"'" for input in inputs]),
            output=" > " + os.path.abspath(output) if output else ""
        )
    )

def cat(
    input,
    output,
    zip=False,
    append=False
    ):

    if not isinstance(input, list):
        inputs=[input]
    else:
        inputs=input

    cat_call = "cat"    
    if zip:
        cat_call = "zcat"

    return Job(
        inputs,
        [output],
        command="""\
{cat} {input} {append}{output}""".format(
            cat=cat_call,
            input=" ".join(inputs) if input else "",
            append=">" if append else "",
            output="> " + output if output else ""
        )
    )

def cut(
    input,
    output,
    options
    ):

    return Job(
        [input],
        [output],
        command="""\
cut {options} {input}{output}""".format(
            options=options,
            input=input if input else "",
            output=" > " + output if output else "",
        )
    )

def paste(
    input,
    output,
    options
    ):

    return Job(
        [input],
        [output],
        command="""\
paste {options} {input}{output}""".format(
            options=options,
            input=input if input else "",
            output=" > " + output if output else "",
        )
    )

def awk(
    input,
    output,
    instructions,
    append=False
    ):

    return Job(
        [input],
        [output],
        command="""\
awk {instructions} {input}{append}{output}""".format(
            instructions=instructions,
            input=input if input else "",
            append=" >" if append else " ",
            output="> " + output if output else "",
        )
    )

def gzip(
    input,
    output,
    ):

    return Job(
        [input],
        [output],
        command="""\
gzip {input}{output}""".format(
            input=input if input else "",
            output=" > " + output if output else "",
        )
    )
