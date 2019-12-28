project('sahap', 'cpp',
  default_options: ['cpp_std=c++11'],
)

inc = include_directories('include')

subdir('src') # libsahap_static_dep
subdir('tests')

executable('sahap', 'main.cpp', dependencies: [
  libsahap_static_dep,
])

