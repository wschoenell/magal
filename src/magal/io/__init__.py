# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

from distutils.version import LooseVersion
from atpy import __version__ as atpy_version
#from magal.io.alhambra import read_ambcat
from magal.io.final_cats import read_set

if LooseVersion(atpy_version) >= LooseVersion('0.9.6'):
    from atpy.registry import register_set_reader #@UnresolvedImport
    from atpy.registry import register_reader #@UnresolvedImport
else:
    from atpy import register_set_reader #@UnresolvedImport @Reimport
    from atpy import register_reader #@UnresolvedImport @Reimport

#register_set_reader('alhambra_catalog', read_ambcat)
register_set_reader('alhambra_catalog', read_set)