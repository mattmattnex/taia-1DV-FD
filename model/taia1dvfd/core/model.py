from dataclasses import dataclass
from .config import Config
from .grid   import Grid
from .state  import State
from .operators import Operators

@dataclass(slots=True)
class Model:
    config: Config
    grid: Grid
    state: State
    operators: Operators