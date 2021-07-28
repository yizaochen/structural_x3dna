from setuptools import setup, find_packages

setup(name='strucpara', 
      version='0.1',
      packages=find_packages(),
      url='https://github.com/yizaochen/structural_x3dna.git',
      author='Yizao Chen',
      author_email='yizaochen@gmail.com',
      license='MIT',
      install_requires=[
          'jupyterlab',
          'pandas',
          'numpy'
      ]
      )
