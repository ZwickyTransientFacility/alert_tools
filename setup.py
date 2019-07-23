from setuptools import setup

setup(name='alert_tools',
      version='0.1',
      description='Get DC mag',
      url='n/a',
      author='M.Gallardo, E.Bellm',
      author_email='mgallardo@astro.ncu.edu.tw, ebellm@uw.edu',
      license='MIT',
      packages=['alert_tools'],
      install_requires=[
          'markdown',
		  'numpy',
		  'pandas',
		  'astropy',
		  'aplpy'
      ],
      zip_safe=False)

