import sys
sys.path.append('/home/tankp/bin/modeller9v8/modlib')
sys.path.append('/home/tankp/bin/modeller9v8/lib/x86_64-intel8')

from modeller import environ, model
env = environ()
mdl = model(env, file = sys.argv[1])
mdl.write_data('PSA', sys.argv[1], 2)
