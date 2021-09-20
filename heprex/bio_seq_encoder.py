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
# Created by dsimm on 14/10/15.
#

# Inclusions
import array
import json
# Load BioPython packages
from Bio.SeqRecord import SeqRecord


# JSON specialized encoder for BioPython BioSeq objects
class BioSeqEncoder(json.JSONEncoder):
    def default(self, obj):
        # Single SeqRecord object
        if isinstance(obj, SeqRecord):
            return {'seq': str(obj.seq), 'id': obj.id, 'name': obj.name, 'description': obj.description,
                    'dbxrefs': obj.dbxrefs}
        # Array of SeqRecord objects
        elif isinstance(obj, array):
            for (i, item) in enumerate(obj):
                json_obj = dict()
                if isinstance(item, SeqRecord):
                    json_obj[i] = {'seq': str(obj.seq), 'id': obj.id, 'name': obj.name, 'description': obj.description,
                        'dbxrefs': obj.dbxrefs}
            return json_obj
        # Let the base class default method raise the TypeError
        return json.JSONEncoder.default(self, obj)
