# **********************************************************************
#  Copyright notice
#
#  (c) 2015-2021 Dominic Simm <dominic.simm@cs.uni-goettingen.de>
#  All rights reserved.
#
#  This file is part of Odysseus.
#
#  Odysseus is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 2 of the License, or
#  (at your option) any later version.
#
#  Odysseus is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with Odysseus.  If not, see <http://www.gnu.org/licenses/>.
# **********************************************************************

# Changes 2015-2021 by Dominic Simm <dominic.simm@cs.uni-goettingen.de>
# See the ChangeLog or git repository for details.

#
# Created by dsimm1 on 04/08/15.
# Updated by dsimm1 on 07/06/19.
#

# Inclusions
from setuptools import setup

setup(name='odysseus',
      version='1.2',
      description='Gene design software enabling selective protein biosynthesis regulation in modern expression systems.',
      url='https://github.com/dsimm/Odysseus',
      author='Dominic Simm',
      author_email='dominic.simm@cs.uni-goettingen.de',
      license='CC BY-NC-SA 4.0',
      packages=['heprex'],
      zip_safe=False)
