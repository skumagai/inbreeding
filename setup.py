try:
    from setuptools import setup
except:
    from distutils.core import setup

config = {
        'description': 'A forward simulator for partially selfing population based on simuPOP',
        'author': 'Seiji Kumagai',
        'url': 'https://github.com/skumagai/selfingsim.git',
        'author_email': 'seiji.kumagai@gmail.com',
        'verson': 1.0,
        'install_requires': ['nose'],
        'packages': ['selfingsim'],
        'scripts': ['bin/selfingsim'],
        'name': 'selfingsim'
        }

setup(**config)
