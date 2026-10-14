from dataclasses import dataclass
from .config import Config
from .grid   import Grid
from .state  import State

@dataclass(slots=True)
class Model:
    config: Config
    grid: Grid
    state: State