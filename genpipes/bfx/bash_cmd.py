################################################################################
# Copyright (C) 2025 C3G, The Victor Phillip Dahdaleh Institute of Genomic Medicine at McGill University
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
    """
    Invokes mkdir command.
    """

    return Job(
        [],
        [folder],
        command=f"mkdir -p {folder} && touch {folder}",
        removable_files=[folder] if remove else []
    )

def chdir(folder):
    """
    Invokes cd command.
    """
    return Job(
        [],
        [folder],
        command=f"cd {folder}"
    )

def ln(
    target_file,
    link,
    input_file=None,
    output=None,
    remove=False
    ):
    """
    Invokes ln command.
    """

    input_files = [input_file] if input_file else [target_file]
    outputs = [output] if output else [link]

    return Job(
        input_files,
        outputs,
        command=f"ln -s -f {target_file} {link}",
        removable_files=[link] if remove else []
    )

def mv(
    source,
    target,
    force=False,
    extra=None
    ):
    """
    Invokes mv command.
    """

    if isinstance(source, list):
        input_files = source
    else:
        input_files = [source]
    return Job(
        input_files,
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
    """
    Invokes cp command.
    """

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
    """
    Invokes rm command.
    """

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
    """
    Invokes touch command.
    """
    return Job(
        [],
        [],
        command=f"touch {target}"
    )

def md5sum(
    input_file,
    output=None,
    check=False
    ):
    """
    Invokes md5sum command.
    """

    if not output:
        output = f"{input_file}.md5"

    return Job(
        [input_file],
        [output],
        command=f"""md5sum {"-c " if check else ""}{input_file} > {output}"""
    )

def cat(
    input_file,
    output,
    zipped=False,
    append=False
    ):
    """
    Invokes cat command.
    """

    if not isinstance(input_file, list):
        input_files=[input_file]
    else:
        input_files=input_file

    cat_call = "cat"
    if zipped:
        cat_call = "zcat"

    return Job(
        input_files,
        [output],
        command="""\
{cat} {input_file} {append}{output}""".format(
            cat=cat_call,
            input_file=" ".join(input_files) if input_file else "",
            append=">" if append else "",
            output="> " + output if output else ""
        )
    )

def cut(
    input_file,
    output,
    options=None
    ):
    """
    Invokes cut command.
    """

    return Job(
        [input_file],
        [output],
        command="""\
cut {options} {input_file}{output}""".format(
            options=options if options else "",
            input_file=input_file if input_file else "",
            output=" > " + output if output else "",
        )
    )

def paste(
    input_file,
    output,
    options=None
    ):
    """
    Invokes paste command.
    """

    return Job(
        [input_file],
        [output],
        command="""\
paste {options} {input_file}{output}""".format(
            options=options if options else "",
            input_file=input_file if input_file else "",
            output=" > " + output if output else "",
        )
    )

def awk(
    input_file,
    output,
    instructions,
    append=False
    ):
    """
    Invokes awk command.
    """

    return Job(
        [input_file],
        [output],
        command="""\
awk {instructions} {input_file}{append}{output}""".format(
            instructions=instructions,
            input_file=input_file if input_file else "",
            append=" >" if append else " ",
            output="> " + output if output else "",
        )
    )

def sed(
    input_file,
    output,
    instructions
    ):
    """
    Invokes sed command.
    """

    return Job(
        [input_file],
        [output],
        command="""\
sed {instructions} {input_file} {output}""".format(
            instructions=instructions,
            input_file=input_file if input_file else "",
            output="> " + output if output else ""
        )
    )

def grep(
    input_file,
    output,
    instructions
    ):
    """
    Invokes grep command.
    """

    return Job(
        [input_file],
        [output],
        command="""\
grep {instructions} {input_file} {output}""".format(
            instructions=instructions,
            input_file=input_file if input_file else "",
            output="> " + output if output else ""
        )
    )

def sort(
    input_file,
    output,
    instructions,
    extra=None
    ):
    """
    Invokes sort command.
    """

    return Job(
        [input_file],
        [output],
        command="""\
sort {instructions} {input_file} {output}{extra}""".format(
            instructions=instructions,
            input_file=input_file if input_file else "",
            output="> " + output if output else "",
            extra=extra if (extra and output) else ""
        )
    )

def zip(
    input_file,
    output,
    recursive=False
    ):
    """
    Invokes zip command.
    """

    input_files = [input_file] if not isinstance(input_file, list) else input_file
    if output:
        outputs = [output]
    else:
        output = f"{input_file}.zip"
        outputs = None
    if recursive:
        rec = "-r"
        if isinstance(input_file, list):
            input_file = os.path.dirname(input_file[0])
    else:
        rec = ""
    return Job(
        input_files,
        outputs,
        command=f"zip {rec} {output} {input_file}"
    )

def gzip(
    input_file,
    output,
    options=None
    ):
    """
    Invokes gzip command.
    """

    return Job(
        [input_file],
        [output],
        command="""\
gzip {options}{input_file}{output}""".format(
            options=options if options else "",
            input_file=input_file if input_file else "",
            output=" > " + output if output else "",
        )
    )

def chmod(file, permission):
    """
    Invokes chmod command.
    """

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
    input_files,
    threads,
    options=None,
    ini_section='pigz'
    ):
    """
    Invokes pigz compression tool.
    """

    return Job(
        input_files=input_files,
        output_files=[s + ".gz" for s in input_files],
        module_entries=[
            [ini_section, 'module_pigz']
        ],
        command="""\
pigz {options} \\
  {nthreads} \\
  {input_files}""".format(
            input_files=" ".join(input_files),
            nthreads=threads,
            options=options if options else ""
        )
    )

def ls(target):
    """
    Invokes ls command.
    """

    return Job(
        [target],
        [],
        command=f"ls {target}"
    )

def tar(
    input_files,
    output,
    options="-c",
    file_list=False,
    input_file_dependency=True,
    ):
    """
    Invokes tar compression tool.

    By default it compresses input_files, but it can be used to extract the content
    of a tar with option "-x" and supplying the tar in the output arg.
        tar {options} -f {output} {input_files}

    Args:
        input_files      list
            Paths/strings of files and folders to compress.
        output      str
            Filename preferably with extension .tar or .tar.gz.
        options     str, default = "-c"
            Arguments.
        file_list   boolean, default = False
            Adds a second output file to the Job object.
            The file is a list of the tar content in plain text.
        input_file_dependency    boolean, defautlt = True
            Adds dependencies to input_file files to the Job object.
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
    if input_file_dependency:
        dependencies = input_files
    else:
        dependencies = [None]
    return Job(
        dependencies,
        outputs,
        command="""\
tar {args} -f {output} {input_files}{stdout}""".format(
            args=args,
            output=output,
            input_files=" ".join(input_files),
            stdout=" ".join([" >", outputs[-1]]) if file_list else ""
        )
    )
