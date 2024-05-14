import shutil
import os
from pathlib import Path

def main():
    builddir = os.path.join('.', '_build', 'html')
    outdir = os.path.join('..', 'docs')

    # If docs dir exists, delete it:
    shutil.rmtree(outdir)

    # Copy over the files from the HTML directory:
    shutil.copytree(builddir, outdir)

    # Remove the `_sources` directory:
    shutil.rmtree(os.path.join(outdir, '_sources'))

    # And include the .nojekyll file:
    Path(os.path.join(outdir, '.nojekyll')).touch()

    print('Docs published to {}.'.format(outdir))


if __name__ == "__main__":
    main()