import glob
import os

def main():
    docpaths = ['../+sonic']
    outdir = './api_src'

    # Generate the stub directory if it doesn't exist:
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    # And process:
    for path in docpaths:
        genStubs(os.path.split(path)[-1], path, outdir)

# Stub generation:
def genStubs(module, path, outdir):
    # Grab all .m-files in the directory of interest:
    mfiles = glob.glob(os.path.join(path, '*.m'))
    
    rst_fn_out = []

    # Generate the individual stubs:
    for file in mfiles:
        classname = os.path.basename(file).split('.')[0]
        rst_fn_out.append(genAutoClassPage(module, classname, outdir))

    # Generate the tocfile, referecing each stubfile:
    template = '''API\n===\n.. toctree::\n'''
    for fn_out in rst_fn_out:
        template = template + '\t{}\n'.format(os.path.join(outdir, fn_out))

    # And write to file:
    with open('api.rst', 'w') as f:
        f.write(template) 


def genAutoClassPage(module, classname, outdir):

    # Template for a basic class dump:
    template = '''{}\n==========================\n\n.. automodule::  {}
    \n\n.. autoclass:: {}\n\t:show-inheritance:\n\t:members:\n\t:undoc-members:\n'''.format(classname, module, classname)

    fn_out = classname.lower()

    # And write to file:
    with open(os.path.join(outdir, fn_out)+'.rst', 'w') as f:
        f.write(template)

    return fn_out

if __name__ == "__main__":
    main()