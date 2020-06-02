"""
AlchemicalITP
The Gromacs itp file parser and writer for multi-step alchemical transformation.
"""

# Add imports here
from .top import Topology

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
