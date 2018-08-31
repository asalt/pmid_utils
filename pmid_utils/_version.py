import os
from datetime import date
from glob import glob

BASE_DIR = os.path.dirname(os.path.realpath(__file__))

def _calculate_latest_edit(files, path=BASE_DIR):
    return max([os.stat(os.path.join(path, f)).st_mtime for f in files])

_package_files = glob('*py')
__version__ = '0.1.040'
__copyright__ = date.fromtimestamp(_calculate_latest_edit(_package_files)).__str__()
