
with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

import re

long_description = re.sub(
    r':py:[a-z]+:`([^`]+)`',
    r'``\1``',
    readme + '\n\n' + history)

long_description = re.sub(
    r':issue:`([^`]+)`',
    r'`GH issue #\1 <https://github.com/ClimateImpactLab/climate_toolbox/issues/\1>`_',
    long_description)

long_description = re.sub(
    r':pull:`([^`]+)`',
    r'`GH PR #\1 <https://github.com/ClimateImpactLab/climate_toolbox/pull/\1>`_',
    long_description)
