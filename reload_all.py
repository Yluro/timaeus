import importlib
import sys
import os
import olx

def _plugin_dir():
    """The folder this file lives in, so reloading works wherever the plugin is
    installed instead of assuming a fixed path under BaseDir(). Olex2 can exec
    this file without setting __file__, so fall back to the path the plugin
    registered when it loaded."""
    try:
        return os.path.dirname(os.path.abspath(__file__))
    except NameError:
        return olx.GetVar('Timaeus_plugin_path')

def reload_all():
    base = _plugin_dir()
    if base not in sys.path:
        sys.path.insert(0, base)

    for f in os.listdir(base):
        if f.endswith('.py') and f != 'reload_all.py':
            module_name = f[:-3]
            if module_name in sys.modules:
                importlib.reload(sys.modules[module_name])
            else:
                importlib.import_module(module_name)

if __name__ == "__main__":
    reload_all()
