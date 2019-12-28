src = [
	'Chromosome.cpp',
	'Genome.cpp',
	'Allele.cpp',
	'InputReader.cpp',
]

eigen_dep = dependency('eigen', fallback: ['eigen', 'eigen_dep'])

libsahap_static = static_library(
	'sahap', src,
	include_directories: inc,
	dependencies: eigen_dep,
)

libsahap_static_dep = declare_dependency(
	link_with: libsahap_static,
	include_directories: inc,
	dependencies: eigen_dep,
)

libsahap_shared = shared_library(
	'sahap', src,
	include_directories: inc,
	dependencies: eigen_dep,
)
