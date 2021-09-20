# **********************************************************************
#  Copyright notice
#
#  (c) 2015-2021 Dominic Simm <dominic.simm@cs.uni-goettingen.de>
#  All rights reserved.
#
#  This file is part of Odysseus.
#
#  Odysseus is free software: you can redistribute it and/or modify
#  it under the terms of the Creative Commons BY-NC-SA 4.0 License as
#  published by the Creative Commons Corporation, either version 2 of the
#  License, or (at your option) any later version.
#
#  Odysseus is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  Creative Commons BY-NC-SA 4.0 License for more details.
#
#  You should have received a copy of the Creative Commons BY-NC-SA 4.0 License
#  along with Odysseus.  If not, see <https://creativecommons.org/licenses>.
# **********************************************************************

# Changes 2015-2021 by Dominic Simm <dominic.simm@cs.uni-goettingen.de>
# See the ChangeLog or git repository for details.

#
# Created by dsimm on 04/08/15.
#

# Load general Python packages
import os.path, pickle, tempfile
from datetime import datetime


def millis(start_time):
    """Helper function to measure the elapsed milliseconds since passed timepoint (e.g. start of the program)

    Args:
        start_time (datetime): timepoint to measure from.

    Returns:
        float: The measured time in milliseconds
    """

    td = datetime.now() - start_time
    return td.total_seconds() * 1000


# Returns current time in passed format
def timestamp():
    """Helper function to

    Args:
        start_time (datetime): timepoint to measure from.

    Returns:
        float: The measured time.
    """

    dt = datetime.now()

    return dt.strftime("%Y-%m-%d %H:%M")


# Returns unique identifier for storing temporary data
def uniq_id():
    """Helper function to create a unique identifier for storing temporary data

    Returns:
        str: The unique string extracted from TemporaryFile module.
    """

    f = tempfile.NamedTemporaryFile(delete=False)
    os.unlink(f.name)

    return os.path.basename(f.name)


# Re-Implementation of deepcopy for all data structures
def deepcopy(obj):
    """Helper function: Primitive Re-Implementation of deepcopy for nested dicts, lists, etc. using the Pickle module

    Args:
        obj: Object to copy.

    Returns:
        obj: The pickle duplicated object.
    """

    return pickle.loads(pickle.dumps(obj))


# Re-Implementation of deepcopy for nested dicts, lists, etc.
# def deepcopy(obj):
#     if isinstance(obj, dict):
#         return {deepcopy(key): deepcopy(value) for key, value in obj.items()}
#     if hasattr(obj, '__iter__'):
#         return type(obj)(deepcopy(item) for item in obj)
#     return obj
