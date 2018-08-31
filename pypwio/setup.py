reqs = [
    'pymatgen>=2018.2.13'
]

setup_args = {
    'name': 'pypwio',
    'version': '0.1.0',
    'description': 'Quantum Espresso pw.x I/O',
    'include_package_data': True,
    'install_requires': reqs,
    'author': 'Austin Zadoks',
    'author_email': 'zooks97@gmail.com',
    'keywords': [
        'materials science',
        'quantum espresso',
        'pw.x'
    ],
    'classifiers': [
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',
        'Programming Language :: Python :: 3.6'
    ]
}