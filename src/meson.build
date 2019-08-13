src = [
	'Chromosome.cpp',
	'Genome.cpp',
	'DNAChar.cpp',
	'Read.cpp',
	'InputReader.cpp',
]

libsahap_static = static_library(
	'sahap', src,
	include_directories: inc,
)

libsahap_static_dep = declare_dependency(
	link_with: libsahap_static,
	include_directories: inc,
)

libsahap_shared = shared_library(
	'sahap', src,
	include_directories: inc,
)
