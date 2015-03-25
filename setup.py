from distutils.core import setup

setup(
    name='boxmodel',
    version='0.2',
    author='Laurent Pouilloux',
    author_email='laurent.pouilloux@inria.fr',
    package_dir={'': 'src/'},
    packages=['boxmodel'],
    scripts=['bin/boxmodel-run'],
    url='https://github.com/lpouillo/boxmodel',
    license='LICENSE.txt',
    description='A python module to compute flux-box models of' + \
        ' isotopic elements in human body.',
    long_description=open('README.txt').read()
)
