target_sources(chargefw_core_objects
        PRIVATE
        atom.cpp
        bond.cpp
        conformer.cpp
        molecule.cpp
        molecule_collection.cpp
        periodic_table.cpp
        position.cpp
        elements.h
        elements.cpp

        PUBLIC
        FILE_SET HEADERS
        BASE_DIRS
        ../../include
        FILES
        ../../include/chargefw/core/atom.h
        ../../include/chargefw/core/bond.h
        ../../include/chargefw/core/conformer.h
        ../../include/chargefw/core/molecule.h
        ../../include/chargefw/core/molecule_collection.h
        ../../include/chargefw/core/position.h
        ../../include/chargefw/core/periodic_table.h
)