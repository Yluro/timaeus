import importlib
import sys
import os
import olx

def reload_all():
    base = os.path.join(olx.BaseDir(), 'util', 'pyUtil', 'PluginLib', 'plugin-SymmetryMeasurements')
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