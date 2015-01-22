from setuptools import setup

setup(description='A forward simulator for partially selfing population based on simuPOP',
      long_description="""
        A forward-in-time population genetic simulator for partially selfing organisms.
        Currently three modes of mating (pure hermaphroditism, androdioecy, and gynodioecy)
        are supported.

        A constant-size population evolves under the infinite alleles or sites model.
        """,
      author='Seiji Kumagai',
      url='https://github.com/skumagai/selfingsim.git',
      author_email='seiji.kumagai@gmail.com',
      version=1.0,
      install_requires=['nose'],
      packages=['selfingsim'],
#      scripts=['scripts/selfingsim'],
      entry_points={
          'console_scripts': [
              'selfingsim=selfingsim:run',
          ],
      },
      name='selfingsim',
      license='GPLv2',
      classifiers=[
         'Development Status :: 5 - Production/Stable',
         'Intended Audience :: Science/Research',
         'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)'
         'Environment :: Console',
         'Programming Language :: Python :: 2.7',
         'Programming Language :: Python :: 3',
         'Topic :: Scientific/Engineering :: Bio-Informatics'
         ],
      keywords='biology simulation population_genetics selfing',
      package_data={
          'example': ['examples/example.json', 'examples/example.json.annotated']
      }
)
