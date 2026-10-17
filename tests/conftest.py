import importlib.util
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parent.parent / "scripts"


def load_script(filename):
    """Import a pipeline script by file name (they start with a step number, so plain import fails)."""
    spec = importlib.util.spec_from_file_location(filename.removesuffix(".py").replace("-", "_"), SCRIPTS / filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module
