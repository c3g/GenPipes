################################################################################
# Copyright (C) 2014, 2023 GenAP, McGill University and Genome Quebec Innovation Centre
#
# This file is part of GenPipes.
#
# GenPipes is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# GenPipes is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with GenPipes.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

# Python Standard Modules
import os

# MUGQIC Modules
from ..core.job import Job

def mkdir(
    folder,
    remove=False
    ):

    return Job(
        [],
        [folder],
        command="""\
mkdir -p {directory} && \\
touch {directory}""".format(
            directory=folder
        ),
        removable_files=[folder] if remove else []
    )

def chdir(folder):
    return Job(
        [],
        [folder],
        command="""\
cd {directory}""".format(
            directory=folder
        )
    )

def ln(
    target_file,
    link,
    input=None,
    output=None,
    remove=False
    ):

    inputs = [input] if input else [target_file]
    outputs = [output] if output else [link]

    return Job(
        inputs,
        outputs,
        command="""\
ln -s -f \\
  {target_file} \\
  {link}""".format(
            target_file=target_file,
            link=link
        ),
        removable_files=[link] if remove else []
    )

def mv(
    source,
    target,
    force=False,
    extra=None
    ):
    if isinstance(source, list):
        inputs = source
    else: 
        inputs = [source]
    return Job(
        inputs,
        [target],
        command="""\
mv {force}{source} \\
   {dest}{extra}""".format(
            force="-f " if force else "",
            source=source,
            dest=f"-t {target}" if os.path.isdir(target) else target,
            extra=extra if extra else ""
        )
    )

def cp(
    source,
    target,
    recursive=False,
    update=False
    ):

    return Job(
        [source],
        [target],
        command="""\
cp {rec}{upd}'{source}' {dest}""".format(
            source=source,
            dest=target,
            rec="-r " if recursive else "",
            upd="-u " if update else ""
        )
    )

def rm(
    source,
    recursive=True,
    force=True
    ):

    return Job(
        [source],
        [],
        command="""\
rm {rec}{force}{source}""".format(
            rec="-r " if recursive else "",
            force="-f " if force else "",
            source=source
        )
    )

def touch(target):
    return Job(
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
    options=None
    ):

    return Job(
        [input],
        [output],
        command="""\
cut {options} {input}{output}""".format(
            options=options if options else "",
            input=input if input else "",
            output=" > " + output if output else "",
        )
    )

def paste(
    input,
    output,
    options=None
    ):

    return Job(
        [input],
        [output],
        command="""\
paste {options} {input}{output}""".format(
            options=options if options else "",
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

def sed(
    input,
    output,
    instructions
    ):

    return Job(
        [input],
        [output],
        command="""\
sed {instructions} {input} {output}""".format(
            instructions=instructions,
            input=input if input else "",
            output="> " + output if output else ""
        )
    )

def grep(
    input,
    output,
    instructions
    ):
    return Job(
        [input],
        [output],
        command="""\
grep {instructions} {input} {output}""".format(
            instructions=instructions,
            input=input if input else "",
            output="> " + output if output else ""
        )
    )

def sort(
    input,
    output,
    instructions,
    extra=None
    ):
    return Job(
        [input],
        [output],
        command="""\
sort {instructions} {input} {output}{extra}""".format(
            instructions=instructions,
            input=input if input else "",
            output="> " + output if output else "",
            extra=extra if (extra and output) else ""
        )
    )

def zip(
    input,
    output,
    recursive=False
    ):
    inputs = [input] if not isinstance(input, list) else input
    if output:
        outputs = [output] 
    else:
        output = f"{input}.zip"
        outputs = None
    if recursive:
        rec = "-r"
        if isinstance(input, list):
            input = os.path.dirname(input[0])
    else:
        rec = ""
    return Job(
        inputs,
        outputs,
        command=f"zip {rec} {output} {input}"
    )

def gzip(
    input,
    output,
    options=None
    ):

    return Job(
        [input],
        [output],
        command="""\
gzip {options}{input}{output}""".format(
            options=options if options else "",
            input=input if input else "",
            output=" > " + output if output else "",
        )
    )

def chmod(file, permission):
    return Job(
        [],
        [file],
        command="""\
chmod {permission} {file}""".format(
            permission=permission,
            file=file
        )
    )

def pigz(
    inputs,
    threads,
    options=None,
    ini_section='pigz'
    ):
    return Job(
        input_files=inputs,
        output_files=[s + ".gz" for s in inputs],
        module_entries=[
            [ini_section, 'module_pigz']
        ],
        command="""\
pigz {options} \\
  {nthreads} \\
  {input_files}""".format(
            input_files=" ".join(inputs),
            nthreads=threads,
            options=options if options else ""
        )
    )

def ls(
    target
    ):
    return Job(
        [target],
        [],
        command="""\
ls {path}""".format(
            path=target
        )
    )

def tar(
    inputs,
    output,
    options="-c",
    file_list=False,
    input_dependency=True,
    ):
    """
    Invokes tar compression tool.

    By default it compresses inputs, but it can be used to extract the content
    of a tar with option "-x" and supplying the tar in the output arg.
        tar {options} -f {output} {inputs}

    Args:
        inputs      list
            Paths/strings of files and folders to compress.
        output      str
            Filename preferably with extension .tar or .tar.gz.
        options     str, default = "-c"
            Arguments.
        file_list   boolean, default = False
            Adds a second output file to the Job object.
            The file is a list of the tar content in plain text.
        input_dependency    boolean, defautlt = True
            Adds dependencies to input files to the Job object.
            Set to False to remove dependency.

    Returns:
        Job object
    """
    outputs = [output]
    if file_list:
        args = " ".join([options, "-vv"])
        outputs.append(
                "".join([os.path.splitext(os.path.splitext(output)[0])[0],
                         ".list"])
                )
    else:
        args = options
    if input_dependency:
        dependencies = inputs
    else:
        dependencies = [None]
    return Job(
        dependencies,
        outputs,
        command="""\
tar {args} -f {output} {inputs}{stdout}""".format(
            args=args,
            output=output,
            inputs=" ".join(inputs),
            stdout=" ".join([" >", outputs[-1]]) if file_list else ""
        )
    )
