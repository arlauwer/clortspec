"""The model configuration, all in this folder:

  config/model.py    the model, shared by all targets (geometry, energies, torus grid)
  config/<name>.py   one target: redshift, Galactic N_H, combined grid (read by load_target)

Scripts do `import config as C`; that gives everything in model.py.
"""

from .model import *  # noqa: F401,F403
