from setuptools import setup, find_packages
install_requires = [
    'setuptools',
    'numpy >= 1.10.0',
    'scipy >= 0.14.0',
    'pynbody >= 0.40',
    'tangos >=1.0.0'
    ]

setup(name='tangos_properties',
      version='1.0.0',
      description='Martin extensions to tangos properties',
      classifiers=[
          "Development Status :: 4 - Beta",
          "Intended Audience :: Developers",
          "Programming Language :: Python :: 3.7",
      ],
      author="Martin Rey",
      author_email="martin.rey.16@ucl.ac.uk",
      license="GNUv3",
      # packages=find_packages(),
      packages=['tangos_properties'],
      entry_points={},
      include_package_data=True,
      zip_safe=False,
      python_requires='>=3.6',
      install_requires=install_requires,
      )